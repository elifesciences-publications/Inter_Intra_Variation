clc;
clear;
oss = struct2dataset(tdfread('../../raw_data/SC_OSSdatatable.txt'));
hp = struct2dataset(tdfread('../../raw_data/SC_HPdatatable.txt'));

HPmids = cellstr(hp.id);
OSSmids = cellstr(oss.mid);
OSScids = cellstr(oss.cid);
% OSSids=strcat('mouse_', cellfun(@(x) x(isstrprop(x,'digit')), cellstr(oss.id), 'UniformOutput', false));

[ossC,osscia,osscic] = unique(OSScids);
[ossmC,ossmia,ossmic] = unique(OSSmids);
[hpmC,hpmia,hpmic] = unique(HPmids);

