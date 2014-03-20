function graphical_user_interface(nnodes,coor,nelem,elnodes,StressNode,...
    U,VonMises,Tresca,nloadinc,All_soln)

function plotelements(coordinates, color, alpha, linestyle)
    if isstr(color)
        c = @(i) color;
    else
        c = @(i) color(i, :);
    end
    for i=1:nelem
        h=fill(coordinates(elnodes(i, 2:5),1),coordinates(elnodes(i, 2:5),2),c(i));
        set(h,'FaceAlpha',alpha);
        set(h, 'LineStyle', linestyle);
    end
end

function [maxstress, minstress] = calcrange(manual, stress, maxstress, minstress)
if ~manual
    maxstress = max(max(stress));
    minstress = min(min(stress));
end
if(maxstress-minstress<2e-6)
    maxstress = maxstress + 1e-6;
    minstress = minstress - 1e-6;
end
caxis([minstress maxstress]);
end

%Graphical output
dpm=[U(1:2:2*nnodes) U(2:2:2*nnodes)];

magfac = 1;

dcoor = coor + magfac.*dpm;

disp('Click on an option in the plot menu window.')
disp('Please be patient in the case of large meshes. It may take a few seconds .......')

choice = 0;
manual = 0;
maxstress = inf;
minstress = 0;
while choice ~= 12
    choice = menu('Plot options','Undeformed shape','Deformed shape', ...
        'Overlay','S11 contour','S22 contour','S12 contour', ...
        'Von Mises contour','Tresca contour', ...
        'Manual colorbar scaling', ...
        'Automatic colorbar scaling','Set magnification factor','Quit');
    if any(choice == [1:8])
        figure(1)
        clf
        hold on
    end
    
    if choice == 1
        plotelements(coor, 'b', 0.8, '-');
        title('Undeformed shape')
    elseif choice == 2
        plotelements(dcoor, 'b', 0.8, '-');
        title(['Deformed shape, magnification factor = ',num2str(magfac)])
    elseif choice == 3
        plotelements(coor, 'b', 0.2, '-');
        if nloadinc==1
            plotelements(dcoor, 'b', 0.8, '-');
        else
            for kk=1:nloadinc
                dpm = All_soln(1+nnodes*(kk-1):nnodes*kk,:);
                dcoor = coor + dpm;
                plotelements(dcoor, 'b', 0.2, '-');
            end
        end
        title('Deformed over undeformed shape')
    elseif choice == 4
        stress = StressNode(:,[2 5 8 11]);
        plotelements(dcoor, StressNode(:, [2, 5, 8, 11]), 1, 'none')
        [maxstress, minstress] = calcrange(manual, stress, maxstress, minstress);
        colorbar('horiz')
        %title('Contour plot: \sigma_{11}')
    elseif choice==5
        stress = StressNode(:, [3, 6, 9, 12]);
        plotelements(dcoor, stress, 1, 'none');
        [maxstress, minstress] = calcrange(manual, stress, maxstress, minstress);
        colorbar('horiz')
        title('Contour plot: \sigma_{22}')
    elseif choice==6
        stress = StressNode(:, [4, 7, 10, 13]);
        plotelements(dcoor, stress, 1, 'none');
        [maxstress, minstress] = calcrange(manual, stress, maxstress, minstress);
        colorbar('horiz')
        title('Contour plot: \sigma_{12}')
    elseif choice==7
        stress = VonMises;
        plotelements(dcoor, stress, 1, 'none');
        [maxstress, minstress] = calcrange(manual, stress, maxstress, minstress);
        colorbar('horiz')
        title('Contour plot: Von Mises')
    elseif choice==8
        stress = Tresca;
        plotelements(dcoor, stress, 1, 'none');
        [maxstress, minstress] = calcrange(manual, stress, maxstress, minstress);
        colorbar('horiz')
        title('Contour plot: Tresca')
    elseif choice == 9
        manual = 1;
        minstress = input('Minimum contour level ? ');
        maxstress = input('Maximum contour level ? ');
        [maxstress, minstress] = calcrange(manual, stress, maxstress, minstress);
        %colorbar('horiz')
        colorbar('vert')
        figure(1)
        refresh
    elseif choice == 10
        manual = 0;
    elseif choice == 11
        magfac = input('Enter new displacement magnification factor. ');
        dcoor = coor + magfac.*dpm;
    end
    if any(choice == 1:8)
        hold off
        axis equal
    end
    
end
end
