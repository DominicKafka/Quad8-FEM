function graphical_user_interface(nnodes,coor,nelem,elnodes,StressNode,...
    U,VonMises,Tresca,nloadinc,All_soln)

%Graphical output
dpm=[U(1:2:2*nnodes) U(2:2:2*nnodes)];
maxdispx = max(dpm(:,1));
mindispx = min(dpm(:,1));
maxdispy = max(dpm(:,2));
mindispy = min(dpm(:,2));
maxcoordx = max(coor(:,1));
mincoordx = min(coor(:,1));
maxcoordy = max(coor(:,2));
mincoordy = min(coor(:,2));

magfacx = 0.1*(maxcoordx-mincoordx)/(maxdispx-mindispx);
magfacy = 0.1*(maxcoordy-mincoordy)/(maxdispy-mindispy);

magfac = 1;

dcoor = coor + magfac.*dpm;

disp('Click on an option in the plot menu window.')
disp(['Please be patient in the case of large meshes. It may take a few seconds .......'])

fac1 = 1;
fac2 = 1;

choice = 0;
manual = 0;
while choice ~= 12
    choice = menu('Plot options','Undeformed shape','Deformed shape', ...
        'Overlay','S11 contour','S22 contour','S12 contour', ...
        'Von Mises contour','Tresca contour', ...
        'Manual colorbar scaling', ...
        'Automatic colorbar scaling','Set magnification factor','Quit');
    if choice == 1
        figure(1)
        clf
        hold on
        for i=1:nelem
            h=patch(coor(elnodes(i,[2 3 4 5 2]),1),coor(elnodes(i,[2 3 4 5 2]),2),'b');
            set(h,'FaceAlpha',0.8)
        end
        hold off
        axis equal
        title('Undeformed shape')
    elseif choice == 2
        figure(1)
        clf
        hold on
        for i=1:nelem
            h=patch(dcoor(elnodes(i,[2 3 4 5 2]),1),dcoor(elnodes(i,[2 3 4 5 2]),2),'b');
            set(h,'FaceAlpha',0.8)
        end
        hold off
        axis equal
        title(['Deformed shape, magnification factor = ',num2str(magfac)])
    elseif choice == 3
        figure(1)
        clf
        hold on
        for i=1:nelem
            h=patch(coor(elnodes(i,[2 3 4 5]),1),coor(elnodes(i,[2 3 4 5]),2),'c');
            set(h,'FaceAlpha',0.2)
        end
        if nloadinc==1
            for i=1:nelem
                h=patch(dcoor(elnodes(i,[2 3 4 5 2]),1),dcoor(elnodes(i,[2 3 4 5 2]),2),'b');
                set(h,'FaceAlpha',0.8)
            end
        else
            for kk=1:nloadinc
                dpm = All_soln(1+nnodes*(kk-1):nnodes*kk,:);
                dcoor = coor + dpm;
                for i=1:nelem
                    h=patch(dcoor(elnodes(i,[2 3 4 5]),1),dcoor(elnodes(i,[2 3 4 5]),2),'b');
                    set(h,'FaceAlpha',0.2)
                end
            end
        end
        hold off
        axis equal
        title('Deformed over undeformed shape')
    elseif choice == 4
        figure(1)
        clf
        hold on
        if manual==0
            maxstress = max(max(StressNode(:,[2 5 8 11])));
            minstress = min(min(StressNode(:,[2 5 8 11])));
        end
        if(maxstress-minstress<2e-6)
            maxstress = maxstress + 1e-6;
            minstress = minstress - 1e-6;
        end
        for i=1:nelem
            x=dcoor(elnodes(i,2:5),1);
            y=dcoor(elnodes(i,2:5),2);
            h=fill(x,y,StressNode(i,[2 5 8 11]));
            set(h,'LineStyle','none')
        end
        hold off
        caxis([minstress maxstress]);
        axis equal
%        colorbar('horiz')
         colorbar('vert')
        %title('Contour plot: \sigma_{11}')
    elseif choice==5
        figure(1)
        clf
        hold on
        if manual == 0
            maxstress = max(max(StressNode(:,[3 6 9 12])));
            minstress = min(min(StressNode(:,[3 6 9 12])));
        end
        if(maxstress-minstress<2e-6)
            maxstress = maxstress + 1e-6;
            minstress = minstress - 1e-6;
        end
        for i=1:nelem
            x=dcoor(elnodes(i,2:5),1);
            y=dcoor(elnodes(i,2:5),2);
            h=fill(x,y,StressNode(i,[3 6 9 12]));
            set(h,'LineStyle','none')
        end
        hold off
        caxis([minstress maxstress]);
        axis equal
        colorbar('horiz')
        title('Contour plot: \sigma_{22}')
    elseif choice==6
        figure(1)
        clf
        hold on
        if manual == 0
            maxstress = max(max(StressNode(:,[4 7 10 13])));
            minstress = min(min(StressNode(:,[4 7 10 13])));
        end
        if(maxstress-minstress<2e-6)
            maxstress = maxstress + 1e-6;
            minstress = minstress - 1e-6;
        end
        for i=1:nelem
            x=dcoor(elnodes(i,2:5),1);
            y=dcoor(elnodes(i,2:5),2);
            h=fill(x,y,StressNode(i,[4 7 10 13]));
            set(h,'LineStyle','none')
        end
        hold off
        caxis([minstress maxstress]);
        axis equal
        colorbar('horiz')
        title('Contour plot: \sigma_{12}')
    elseif choice==7
        figure(1)
        clf
        hold on
        if manual == 0
            maxstress = max(max(VonMises(:,:)));
            minstress = min(min(VonMises(:,:)));
        end
        if(maxstress-minstress<2e-6)
            maxstress = maxstress + 1e-6;
            minstress = minstress - 1e-6;
        end
        for i=1:nelem
            x=dcoor(elnodes(i,2:5),1);
            y=dcoor(elnodes(i,2:5),2);
            h=fill(x,y,VonMises(i,:));
            set(h,'LineStyle','none')
        end
        hold off
        caxis([minstress maxstress]);
        axis equal
        colorbar('horiz')
        title('Contour plot: Von Mises')
    elseif choice==8
        figure(1)
        clf
        hold on
        if manual == 0
            maxstress = max(max(Tresca(:,:)));
            minstress = min(min(Tresca(:,:)));
        end
        if(maxstress-minstress<2e-6)
            maxstress = maxstress + 1e-6;
            minstress = minstress - 1e-6;
        end
        for i=1:nelem
            x=dcoor(elnodes(i,2:5),1);
            y=dcoor(elnodes(i,2:5),2);
            h=fill(x,y,Tresca(i,:));
            set(h,'LineStyle','none')
        end
        hold off
        caxis([minstress maxstress]);
        axis equal
        colorbar('horiz')
        title('Contour plot: Tresca')
    elseif choice == 9
        oldmin = minstress;
        oldmax = maxstress;
        manual = 1;
        minstress = input('Minimum contour level ? ');
        maxstress = input('Maximum contour level ? ');
        if(maxstress-minstress<2e-6)
            maxstress = maxstress + 1e-6;
            minstress = minstress - 1e-6;
        end
        caxis([minstress maxstress]);
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
end
