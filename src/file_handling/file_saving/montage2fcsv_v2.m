function montage2fcsv_v2(MONTAGE_SOURCE, save_file,input_orientation,atlas)
% MONTAGE_SOURCE path to electrodes.xls (xlsx, csv)
%                or cell with header {name1 R A S; name2 R A S; ...}

% save_file fullpath to fcsv file (c:\electrodes.fcsv)
% input_orientation 'ras'
% atlas... text notice of electrodes corresponding to electrodes
% {'abc';'xyz';...} size(atlas,1) == size(MONTAGE_SOURCE,1) -1

if ischar(MONTAGE_SOURCE)
    [~,~,C]=xlsread(MONTAGE_SOURCE);
elseif iscell(MONTAGE_SOURCE)
    C=MONTAGE_SOURCE;
end



TAB=cell(size(C,1),14);
for i=1:size(C,1)
    if isempty(atlas)
        TAB(i,:)=[{['vtkMRMLMarkupsFiducialNode_' num2str(i)]},C(i,2:4),num2cell([0 0 0 1 1 0 1]),C{i,1} {''} {''}];
    else
        if isempty(atlas{i}); atlas{i}='out of HEAD'; end
        TAB(i,:)=[{['vtkMRMLMarkupsFiducialNode_' num2str(i)]},C(i,2:4),num2cell([0 0 0 1 1 0 1]),C{i,1} atlas{i} {''}];
    end
end


if nargin<3
    input_orientation='ras';
end

tidx=sort([findstr(save_file,'/') findstr(save_file,'\')]);
if ~isempty(tidx)
    save_path=save_file(1:tidx(end));
    [~,~,~]=mkdir(save_path);
end

fid=fopen(save_file,'w+');
if strcmpi(input_orientation,'ras')
        fprintf(fid,'# Markups fiducial file version = 4.5\n# CoordinateSystem = 0 \n# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');
elseif strcmpi(input_orientation,'las')
        fprintf(fid,'# Markups fiducial file version = 4.5\n# CoordinateSystem = 1 \n# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');
end

for i=1:size(TAB,1)
    for j=1:size(TAB,2)
        val=cell2mat(TAB(i,j));
        if isnumeric(val)
           if val==0 || val==1
               val=num2str(val);
           else
               val=num2str(val,'%.2f');
           end
        end
        fprintf(fid,[val,',']);
    end
    fprintf(fid,'\n');
end

fclose(fid);
disp('--- fiducial list was created:')
disp(save_file)

