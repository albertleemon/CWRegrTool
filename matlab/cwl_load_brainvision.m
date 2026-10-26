function EEG = cwl_load_brainvision(vhdr)
% CWL_LOAD_BRAINVISION  Read a BrainVision recording (.vhdr/.vmrk/.eeg) into an
% EEGLAB-style EEG structure, without EEGLAB. Same result as bva-io's pop_loadbv for
% continuous binary data: data in single precision scaled by the channel resolution
% (µV; the channel unit field is honoured), events with type = marker description (e.g. 'R128'), 'New Segment' and
% 'DC Correction' as 'boundary', latency in samples (1-based).
%
%   EEG = cwl_load_brainvision('D:\data\sub-01_task-rest_eeg.vhdr')
%
% If the header names a data or marker file that does not exist (e.g. the recording
% was renamed after export), the .eeg/.vmrk with the same name as the header is used.

[vdir, stem] = fileparts(vhdr);
hdr = read_ini(vhdr);
ci  = hdr.commoninfos;
assert(strcmpi(ci.dataformat, 'binary'), 'only binary BrainVision data is supported');
datafile = resolve(vdir, ci.datafile, [stem '.eeg'], 'data');
switch lower(hdr.binaryinfos.binaryformat)
    case 'int_16',        fmt = 'int16';   bps = 2;
    case 'uint_16',       fmt = 'uint16';  bps = 2;
    case 'ieee_float_32', fmt = 'float32'; bps = 4;
    otherwise, error('unsupported binary format %s', hdr.binaryinfos.binaryformat);
end
nch = str2double(ci.numberofchannels);

EEG = empty_eeg();
EEG.comments = ['Original file: ' ci.datafile];
EEG.srate    = 1e6 / str2double(ci.samplinginterval);
EEG.nbchan   = nch;
scale = ones(1, nch);
z = repmat({[]}, 1, nch);
EEG.chanlocs = struct('labels', z, 'ref', '', 'theta', z, 'radius', z, 'X', z, 'Y', z, 'Z', z, ...
                      'sph_theta', z, 'sph_phi', z, 'sph_radius', z, 'type', '');
for c = 1:nch
    f = strsplit(hdr.channelinfos{c}, ',', 'CollapseDelimiters', false);
    EEG.chanlocs(c).labels = strrep(f{1}, '\1', ',');
    if numel(f) >= 3 && ~isnan(str2double(f{3})), scale(c) = str2double(f{3}); end
    if numel(f) >= 4, scale(c) = scale(c) * unit_to_uv(f{4}, EEG.chanlocs(c).labels); end
    % [Coordinates] radius,theta,phi (BrainVision) -> EEGLAB spherical / Cartesian / polar,
    % the conversions of pop_loadbv + pop_chanedit
    if isfield(hdr, 'coordinates') && c <= numel(hdr.coordinates) && ~isempty(hdr.coordinates{c})
        v = str2double(strsplit(hdr.coordinates{c}, ','));
        if numel(v) == 3 && any(v ~= 0)
            st = v(3) - 90 * sign(v(2));  sp = -abs(v(2)) + 90;  r = v(1);
            [x, y, z] = sph2cart(st / 180 * pi, sp / 180 * pi, r);
            EEG.chanlocs(c).sph_theta = st;  EEG.chanlocs(c).sph_phi = sp;  EEG.chanlocs(c).sph_radius = r;
            EEG.chanlocs(c).X = x;  EEG.chanlocs(c).Y = y;  EEG.chanlocs(c).Z = z;
            EEG.chanlocs(c).theta = -st;  EEG.chanlocs(c).radius = 0.5 - sp / 180;
        end
    end
end

fid = fopen(datafile, 'r', 'ieee-le');
assert(fid > 0, 'cannot open %s', datafile);
fseek(fid, 0, 'eof');  npts = floor(ftell(fid) / (nch * bps));  fseek(fid, 0, 'bof');
switch lower(ci.dataorientation)
    case 'multiplexed', EEG.data = fread(fid, [nch, npts], [fmt '=>float32']);
    case 'vectorized',  EEG.data = fread(fid, [npts, nch], [fmt '=>float32']).';
    otherwise, error('unsupported data orientation %s', ci.dataorientation);
end
fclose(fid);
for c = 1:nch
    EEG.data(c, :) = EEG.data(c, :) * scale(c);
end
EEG.pnts  = size(EEG.data, 2);
EEG.trials = 1;
EEG.xmin  = 0;
EEG.xmax  = (EEG.pnts - 1) / EEG.srate;
EEG.times = (0:EEG.pnts - 1) / EEG.srate * 1000;
EEG.ref   = 'common';

if isfield(ci, 'markerfile')
    mrk = read_ini(resolve(vdir, ci.markerfile, [stem '.vmrk'], 'marker'));
    keys = fieldnames(mrk.markerinfos_raw);
    ev = struct('latency', {}, 'duration', {}, 'channel', {}, 'bvtime', {}, 'bvmknum', {}, ...
                'visible', {}, 'type', {}, 'code', {}, 'urevent', {});
    for i = 1:numel(keys)
        f = strsplit(mrk.markerinfos_raw.(keys{i}), ',', 'CollapseDelimiters', false);
        f(end+1:6) = {''};
        e.latency  = str2double(f{3});
        e.duration = str2double(f{4});
        e.channel  = str2double(f{5});
        e.bvtime   = str2double(f{6});  if isnan(e.bvtime), e.bvtime = []; end
        e.bvmknum  = str2double(keys{i}(3:end));
        e.visible  = '';
        if any(strcmpi(f{1}, {'New Segment', 'DC Correction'}))
            e.type = 'boundary';  e.duration = NaN;
        else
            e.type = strrep(f{2}, '\1', ',');
        end
        e.code    = f{1};
        e.urevent = i;
        ev(end+1) = e; %#ok<AGROW>
    end
    ev = ev([ev.latency] >= 1 & [ev.latency] <= EEG.pnts);
    EEG.event   = ev;
    EEG.urevent = rmfield(ev, 'urevent');
end
end

function k = unit_to_uv(u, label)
% factor from the channel's unit (4th field of [Channel Infos]) to microvolts; an empty
% unit means µV (BrainVision default). Anything spelled with a micro sign / "u" is µV.
u = strtrim(u);  k = 1;
switch lower(u)
    case 'v',  k = 1e6;
    case 'mv', k = 1e3;
    case 'nv', k = 1e-3;
    otherwise
        if ~isempty(u) && ~strcmpi(u(end), 'v')
            warning('cwl:unit', 'channel %s: unknown unit "%s", treated as microvolts', label, u);
        end
end
end

function f = resolve(vdir, name, fallback, what)
f = fullfile(vdir, strtrim(name));
if ~exist(f, 'file')
    alt = fullfile(vdir, fallback);
    assert(exist(alt, 'file') == 2, 'header names %s file %s, which is missing (no %s either)', what, name, fallback);
    fprintf('header names %s file %s, which does not exist; using %s\n', what, name, fallback);
    f = alt;
end
end

function S = read_ini(file)
% sections -> lower-case field names without spaces; [Channel Infos] -> cell by Ch number,
% [Marker Infos] -> markerinfos_raw.MkN
S = struct();  sec = '';
lines = regexp(fileread(file), '\r?\n', 'split');
for i = 1:numel(lines)
    l = strtrim(lines{i});
    if isempty(l) || l(1) == ';', continue; end
    if l(1) == '['
        sec = lower(regexprep(l(2:end-1), '\s', ''));  continue;
    end
    k = find(l == '=', 1);
    if isempty(k) || isempty(sec), continue; end
    key = strtrim(l(1:k-1));  val = l(k+1:end);
    switch sec
        case {'channelinfos', 'coordinates'}
            S.(sec){str2double(key(3:end))} = val;
        case 'markerinfos'
            S.markerinfos_raw.(key) = val;
        case {'commoninfos', 'binaryinfos'}
            S.(sec).(lower(regexprep(key, '[ :]', '_'))) = val;
    end
end
end

function EEG = empty_eeg()
% field layout of EEGLAB's eeg_emptyset
f = {'setname','filename','filepath','subject','group','condition','session','comments', ...
     'nbchan','trials','pnts','srate','xmin','xmax','times','data','icaact','icawinv', ...
     'icasphere','icaweights','icachansind','chanlocs','urchanlocs','chaninfo','ref','event', ...
     'urevent','eventdescription','epoch','epochdescription','reject','stats','specdata', ...
     'specicaact','splinefile','icasplinefile','dipfit','history','saved','etc','datfile','run','roi'};
EEG = cell2struct(repmat({[]}, numel(f), 1), f, 1);
EEG.setname = 'cwl_pipeline';  EEG.saved = 'no';  EEG.etc = struct();
EEG.event = [];  EEG.urevent = [];
EEG.chaninfo = struct('plotrad', [], 'shrink', [], 'nosedir', '+X', 'icachansind', [], 'removedchans', []);
EEG.reject = struct('rejmanual', []);  EEG.stats = struct('jp', []);
end
