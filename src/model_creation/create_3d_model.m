%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monadic operations on MRI slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(''); disp('Performing monadic operations on MRI slices... Please wait');

monadic_operations

disp('Finished performing monadic operations.')


%% Creation of the vertex model
p = patch(isosurface(mri_t1));
[faces,verts] = reducepatch(p, 0.05);


fig = fig + 1; figure(fig); clf; 
isosurface(mri_t1); 
xlabel('I R/L'); 
ylabel('J A/P'); 
zlabel('K S/I')

% Transform MRI with the T_matrix
ones_column = ones(size(verts,1),1);

% ======================
% funkce isosurface přehazuje pořadí dimenzí, IJK->JIK, proto je potřeba
% zaměnit 1. a 2. dimenzi
Tisosurface=[0 1 0 0;
             1 0 0 0;
             0 0 1 0;
             0 0 0 1];
% STL soubor neobsahuje informaci o orientaci dimezí. 3D slicer proto
% automaticky předpokládá, že STL je v LPS orientaci.
Tras2lps=[-1 0 0 0;
          0 -1 0 0;
          0 0  1 0;
          0 0 0  1];

% tedy je potřeba vertexy převást v pořadí
% 1. Tisousurvace zamění dimenze z JIK->IJK
% 2. T_matrix z MRI udělá přepočet z IJK->RAS
% 3. Tras2lps převede vertexy RAS->LPS
% řadí se od zadu: T3 * T2 * T1 * vertexy

verts = Tras2lps * T_matrix * Tisosurface * [verts ones_column]';
% ======================
verts = verts(1:3,:)';

V_tri = triangulation(faces, verts);

fig = fig + 1; figure(fig); clf; 
trimesh(V_tri); 
axis image; 
colormap gray; 
title("Mesh of the final head model.")

% Prompt the user for a file name and location
save_stlfile_user(V_tri, '3D Model');
