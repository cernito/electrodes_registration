function [ELECTRODE,SURFACE,fh]=elc_read(elc_path)

% read ELC file with HD-EEG contact positions and check overlaped
% electrodes
%
% INPUTS:
%    elc_path... string path to elc-file
%
% OUTPUTS:
%    ELECTRODE.units ... string
%    ELECTRODE.pos... [n x 3], cordinates in space, where n is nuber contacts
%    ELECTRODE.labels... {n x 1}, string cell of contact name
%
%    SURFACE.units ... string
%    SURFACE.pos... [n x 3], cordinates in space, where n is nuber of
%                  surface points
%
%    fh ... figure hendel. If nargout<3 and elc-file is corect, figure is hide



% clc; close all; clear all;
% 
% pacient='PS81/H81';
% cesta='D:/EEG/DATA/';
% 
% 
% elc_file=dir(fullfile(cesta,pacient,'montage','*.elc'));
% fo1=fopen(fullfile(cesta,pacient,'montage',elc_file.name));


% načítání z z ascii souboru (.elc, .txt)
fo1=fopen(elc_path);

if fo1==-1
   error('File not found')
end

% ---------------------------------------------------------------------------------------
% Read data from elc-file
% ---------------------------------------------------------------------------------------
HDR=textscan(fo1,'%s%f',1); % "num of points" <-- N: this is the first line in file
n_electrodes=HDR{1,2};

HDR=textscan(fo1,'%s%s',1); % "units" <-- (mm) second line
ELECTRODE.units=HDR{1,2}{1};

HDR=textscan(fo1,'%s',1); % "position" <-- third line, next N lines after this are electrode labels and positions
ELECTRODE.pos = textscan(fo1,'%*s%n%n%n',n_electrodes,'HeaderLines',0,'Delimiter',{':','\t','\t','\t'},'MultipleDelimsAsOne',1);
ELECTRODE.pos=cell2mat(ELECTRODE.pos); % <-- all x, y, z positions of electrodes in matrix

HDR=textscan(fo1,'%s',1);
ELECTRODE.labels = textscan(fo1,'%s',n_electrodes,'Delimiter',{'\t'},'MultipleDelimsAsOne',1); % reads all labelsinto 'cell array'
ELECTRODE.labels=ELECTRODE.labels{1}; % <-- all electore labels in 'cell array'

HDR=textscan(fo1,'%s%f',1); % num of head points <-- M
n_scalp=HDR{1,2};

HDR=textscan(fo1,'%s%s',1); % units
SURFACE.units=HDR{1,2}{1};

HDR=textscan(fo1,'%s',1); % position
SURFACE.pos = textscan(fo1,'%n%n%n',n_scalp,'HeaderLines',0,'Delimiter',{':','\t','\t','\t'},'MultipleDelimsAsOne',1);
SURFACE.pos=cell2mat(SURFACE.pos); % <--- next M lines are scalp positions, x, y, z?

% ---------------------------------------------------------------------------------------
% Výpočet prùměrné vzdálenosti mezi elktrodama
% ---------------------------------------------------------------------------------------
avr_dist=zeros(size(ELECTRODE.pos,1),1);
for i=1:size(ELECTRODE.pos,1)
    DIST=sqrt(sum((repmat(ELECTRODE.pos(i,:),size(ELECTRODE.pos,1),1)-ELECTRODE.pos).^2,2));
    dist_sort=sort(DIST(DIST>0));
    avr_dist(i,1)=median(dist_sort(1:3)); % úhlopříčka pravoúhlého rastru + 20%
end
avr_dist=mean(avr_dist);

disp(['average distance: ' num2str(avr_dist,'%.2f') ' mm'])

% ---------------------------------------------------------------------------------------
% testuje, zdali nejsou dvě elktrody blízko sebe <0.5 avr_dist 
% ---------------------------------------------------------------------------------------
close_el=[];
for i=1:size(ELECTRODE.pos,1)
    DIST=sqrt(sum((repmat(ELECTRODE.pos(i,:),size(ELECTRODE.pos,1),1)-ELECTRODE.pos).^2,2));
    idx=find(DIST<0.5*avr_dist);
    if length(idx)>1
        close_el=[close_el; sort(idx(:)')];
    end
end

% ---------------------------------------------------------------------------------------
% kontrola překryvu elektrod, nalezení konkrétních názvů
% ---------------------------------------------------------------------------------------
if size(close_el,1)>0
    
    C=corr([close_el'; zeros(3,size(close_el,1))]); C=tril(C.*(1-eye(size(C,1))));
    [ro,co]=find(C==1);
    
    
    disp(['average electrode distance: ' num2str(avr_dist,'%.1f') 'mm'])
    el_text=[];
    for j=1:length(ro)
        dist=sqrt(sum((ELECTRODE.pos(close_el(ro(j),1),:)-ELECTRODE.pos(close_el(ro(j),2),:)).^2));
        el_text=[el_text, ELECTRODE.labels{close_el(ro(j),1)}, '-' ELECTRODE.labels{close_el(ro(j),2)}, ' (' num2str(dist,'%.1f') 'mm)' ', '];
    end
    el_text(end-1:end)=[];
    warning(['ELECTRODE: ' el_text])
    
    fh=figure();
    plot3(ELECTRODE.pos(:,1),ELECTRODE.pos(:,2),ELECTRODE.pos(:,3),'.'); axis image
    title(['ERROR: ' el_text])
    hold on
    for i=1:size(ELECTRODE.labels,1)
        if ~isempty(intersect(i,close_el(ro,:)))
            text(ELECTRODE.pos(i,1),ELECTRODE.pos(i,2),ELECTRODE.pos(i,3),ELECTRODE.labels{i,1},'EdgeColor',[1 0 0])
        else
            text(ELECTRODE.pos(i,1),ELECTRODE.pos(i,2),ELECTRODE.pos(i,3),ELECTRODE.labels{i,1})
        end
    end
end

% ---------------------------------------------------------------------------------------
% zobrazení elektrod v 3D plotu
% ---------------------------------------------------------------------------------------
if nargout>2 && exist('fh','var')==0
    fh=figure();
    plot3(ELECTRODE.pos(:,1),ELECTRODE.pos(:,2),ELECTRODE.pos(:,3),'.'); axis image
    hold on
    for i=1:size(ELECTRODE.labels,1)
        text(ELECTRODE.pos(i,1),ELECTRODE.pos(i,2),ELECTRODE.pos(i,3),ELECTRODE.labels{i,1})
    end
end

fclose(fo1);
