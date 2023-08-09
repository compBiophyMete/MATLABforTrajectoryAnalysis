% The first line of a pdb file is comprised of a header, which needs to
% be removed for the calculation.
pdb = fopen('protein.pdb', 'r');
fgetl(pdb);
buffer = fread(pdb, Inf);
fclose(pdb);
pdb = fopen('fullmodelmin_woheaders.pdb', 'w');
fwrite(pdb, buffer);
fclose(pdb);

%Read the pdb file that does not include any headers!
pdb_woheader = fopen('fullmodelmin_woheaders.pdb', 'r');

%Please, check out the following website https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
pdbFormat='%4s  %5d %4s%4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s';
A=textscan(pdb_woheader,pdbFormat);
x=[A{1,7}];
y=[A{1,8}];
z=[A{1,9}];
atoms = [A{1,3}];
if length(x) == length(y) && length(y) == length(z)
    fprintf('The coordinates has been succesfully read!\n')
    totalNumberOfAtoms = length(x)-1;
else
    fprintf('A severe problem has raisen!\n')
end

% Determining the indeces of the Calpha atoms.
Calpha_atoms_windex = [];
for i=1:length(atoms)
    if atoms(i) == "CA"
        Calpha_atoms_windex = [Calpha_atoms_windex, i];
    else
        continue
    end
end
Calpha_atoms = Calpha_atoms_windex';

%Generating a matrix that includes the coordinates of the Calpha atoms.
coordinates_of_Calpha_atoms = [];
i = 1;
while i <= length(Calpha_atoms)
    coordinates_of_Calpha_atoms = [coordinates_of_Calpha_atoms; x(Calpha_atoms(i)), y(Calpha_atoms(i)), z(Calpha_atoms(i))];
    i = i + 1;
end

%Measuring the pairwaise distances among the Calpha atoms within the read PDB file.
i = 1;
distance_matrix = [];
for i=1:length(coordinates_of_Calpha_atoms)
    for j=1:length(coordinates_of_Calpha_atoms)
    distance_matrix = [distance_matrix;sqrt(power((coordinates_of_Calpha_atoms(i,1) - coordinates_of_Calpha_atoms(j,1)),2) + power((coordinates_of_Calpha_atoms(i,2)- coordinates_of_Calpha_atoms(j,2)),2) + power((coordinates_of_Calpha_atoms(i,3)- coordinates_of_Calpha_atoms(j,3)),2))];
    end
end

%Making the distance matrix n by n matrix
nbyn_distance_matrix = reshape(distance_matrix,[length(Calpha_atoms), length(Calpha_atoms)]);

%Illustrating the distances between atom-pairs
h_distance = heatmap(nbyn_distance_matrix);
xlabel('Residue Number')
ylabel('Residue Number')
XYLabels = 1:length(Calpha_atoms);
string_XYLabels = string(XYLabels);
string_XYLabels(mod(XYLabels,10) ~= 0 ) = " ";
h_distance.XDisplayLabels = string_XYLabels;
h_distance.YDisplayLabels = string_XYLabels;

% Reading the dcd file by using a function embedded in the MDToolBox package
trajectory = readdcd('protein.dcd');
sizeOftrajectory = size(trajectory);
numberOfcoordinates = sizeOftrajectory(2);
numberOfframes = sizeOftrajectory(1);
%X:3n-2, Y:3n-1, Z:3n where n = 1:1:3134
x_trajectory = []; y_trajectory = [] ; z_trajectory = []; organized_trajectory = [];
nframes = 1;
while nframes <= numberOfframes
    for n=1:1:numberOfcoordinates/3
        x_trajectory = [x_trajectory;trajectory(nframes,3*n-2)];
        y_trajectory = [y_trajectory;trajectory(nframes,3*n-1)];
        z_trajectory = [z_trajectory;trajectory(nframes,3*n)];
    end
    fprintf('Progress:%f\n', 100*nframes/numberOfframes)
    nframes = nframes+1;
end
organized_trajectory = [organized_trajectory;x_trajectory, y_trajectory, z_trajectory];
% Arranging the Calpha indeces by referring to the
% obtained indeces in the above
Calpha_atoms_trajindex = []; Calpha_traj = [];
for k = 0:numberOfframes-1
    [Calpha_atoms_trajindex] = [Calpha_atoms_trajindex; Calpha_atoms + k*totalNumberOfAtoms];
end
for i = 1:length(Calpha_atoms_trajindex)
    [Calpha_traj] = [Calpha_traj;organized_trajectory(Calpha_atoms_trajindex(i),1),organized_trajectory(Calpha_atoms_trajindex(i),2), organized_trajectory(Calpha_atoms_trajindex(i),3)];
end
%Inserting a column that shows the residue number
row_res_number=[];
for k = 1:numberOfframes
    row_res_number = [row_res_number,1:length(Calpha_atoms)*k/k];
end
column_res_number = (row_res_number)';
Calpha_traj_wresnumber = [Calpha_traj column_res_number];

flag = -1;
while flag < 0
    % In order for user to start search, one of the specified entries have
    % to given as an input.
    command = input('If you would like to start/contiue your search, please enter sc!. If not, please enter t\n', 's');
    if command == 'sc'
        % If an user enters an invalid atom number, that user will be
        % warned!
     atomnumber1 = input('Enter the first atom:\n');
     atomnumber2 = input('Enter the second atom:\n');
     if atomnumber1 > length(Calpha_atoms) || atomnumber2 > length(Calpha_atoms)
         disp('Invalid atom number')
         break
     end
     % Users have to enter an integer!
     if isnan(atomnumber1) || fix(atomnumber1) ~= atomnumber1 && isnan(atomnumer2) || fix(atomnumber2) ~= atomnumber2
        disp('Please enter an integer')
     end
     distanceOfInterest = [];
     coordinates_atomnumber1=[];
     coordinates_atomnumber2=[];
     for i = 1 : length(Calpha_traj_wresnumber)
        if Calpha_traj_wresnumber(i,4) == atomnumber1
            coordinates_atomnumber1 = [coordinates_atomnumber1; Calpha_traj_wresnumber(i,:)];
        end
        Calpha_traj_wresnumber(i,4) == atomnumber2;
        coordinates_atomnumber2 = [coordinates_atomnumber2; Calpha_traj_wresnumber(i,:)];
     end
     i = 1;
     uperdeterminants = size(coordinates_atomnumber1);
     upperboundary = uperdeterminants(1);
     while i <= upperboundary
         distanceOfInterest = [distanceOfInterest, sqrt(power((coordinates_atomnumber1(i,1) - coordinates_atomnumber2(i,1)),2) + power((coordinates_atomnumber1(i,2)- coordinates_atomnumber2(i,2)),2) + power((coordinates_atomnumber1(i,3)- coordinates_atomnumber2(i,3)),2))];
         i = i +1;
    end
    %Histogram
    %The cut-off has to be numerıcal
    cutoff = input('Please enter the cut-off value:\n');
     if isnan(cutoff) || fix(cutoff) ~= cutoff 
     else
         fprintf('Please, enter a numerical value');
         continue
     end
     %If the standard deviation of the measured distance is above the
     %cut-off, generate the histogram plot.
     if std(distanceOfInterest) > cutoff
        figure(2);histogram(distanceOfInterest, 'Normalization', 'probability');
        xlabel('Distance (Å)')
        ylabel('Probability')
        title(['The distance between Calpha no: ',num2str(atomnumber1), ' and ', num2str(atomnumber2), ', where the cut-off=', num2str(cutoff)])
     else 
         %If the standard devitation of the measure distance is below the
         %cut-off, do not generate the histogram plot!
         continue
     end
     % If a user wants to terminate the search, t has to be entered!
    elseif command == 't'
        display('The search is terminated')
        break
        % In the case of entering the wrong commands, the search gets
        % terminated.
    elseif command ~= 't' || command ~= 'sc'
        display('Please, enter a suitable command')
        break
    end
end


    

      
     





        
    
