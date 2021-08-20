%% Bayesian Inversion Stress Microscopy
%% MATLAB script
%% Vincent NIER, 15/11/16
%% Philippe MARCQ, 08/06/18

%If you use BISM please cite the following reference in your work: 
%Nier, V. et al., Inference of internal stress in a cell monolayer, Biophysical Journal, 110(7), 1625-1635, 2016

% Input: Traction_field.mat, traction force data
% Output: Stress_BISM.mat, internal stress tensor
% Parameters: C, number of data points along the x direction, between xmin and xmax
%             R, number of data points along the y direction, between ymin and ymax		      
%		     kM, number of time steps				
% The regularization parameter Lambda can be estimated by either (meth_Lambda)
%             Maximal A Posteriori (MAP) optimization
%             L-curve method 
% Free stress boundary conditions can be imposed in the prior (BC).

% coarse-graining over domains of fixed size

clear all;
close all;

%% Traction force data
nassay = 6;
ForceName=['Traction_field_control_' num2str(nassay, '%d') '.mat'];
kM=4; % number of time steps
% ForceName='Before_Caly_5.mat';
% kM=2; % number of time steps
Rtrue=58; % number of rows (y direction)
Ctrue=59; % number of columns (x direction)

% %% define central region for spatial averaging
Cmintrue = 20;
Cmaxtrue = 40;
Rmintrue = 20;
Rmaxtrue = 40;


%% Define spatial grid (rectangular)
coeff=5.17392; % conversion coefficient from pixels to micrometers
xmin=coeff*1;
xmax=coeff*Ctrue;
ymin=coeff*1;
ymax=coeff*Rtrue;
xmincenter=coeff*Cmintrue;
xmaxcenter=coeff*Cmaxtrue;
ymincenter=coeff*Rmintrue;
ymaxcenter=coeff*Rmaxtrue;

deltax=(xmax-xmin)/(Ctrue-1); %uniform step of the rectangular grid on x
deltay=(ymax-ymin)/(Rtrue-1); %uniform step of the rectangular grid ox y
xcenter = linspace(xmincenter+deltax, xmaxcenter-deltax, Cmaxtrue-Cmintrue+1); %x coordinate (column from left to right)
ycenter = linspace(ymincenter+deltay, ymaxcenter-deltay, Rmaxtrue-Rmintrue+1);



% coarse-graining scale (must be odd)
nR = 1;
nC = 1;
% number of grid points between two points of the new grid
nR2 = (nR-1)/2;
nC2 = (nC-1)/2;

% define central region after coarse-graining
if nR2 == 0   
    R = Rtrue;    
    Rmin = Rmintrue;
    Rmax = Rmaxtrue;
else
    R = floor(Rtrue/nR2)-1;    
    Rmin = floor(Rmintrue/nR2)-1;
    Rmax = floor(Rmaxtrue/nR2)-1;
end

if nC2 == 0
    C = Ctrue;    
    Cmin = Cmintrue;
    Cmax = Cmaxtrue;
else
    C = floor(Ctrue/nC2)-1;    
    Cmin = floor(Cmintrue/nC2)-1;
    Cmax = floor(Cmaxtrue/nC2)-1;
end


% arrays for central pressure and central shear stress
pcenter = zeros(1,kM);
scenter = zeros(1,kM);

% system size for data analysis
N=R*C;    
% system size for stress tensor field
Ninf=4*N+2*(C+R); 
 
% spatial resolution 
lx=(xmax-xmin)/(C-1); 
ly=(ymax-ymin)/(R-1); 
l = (lx+ly)/2.;

x=xmin:lx:xmax; % x coordinate (column from left to right)
y=ymin:ly:ymax; % y coordinate (row from top to bottom)
 
%% Parameters
BC=0; 
% 1: free boundary conditions \sigma_{ik}n_k=0
% 0: otherwise
meth_Lambda=3; 
% method to determine the regularization parameter Lambda
% 1: MAP 
% 2: L-curve 
% 3: fixed value
noise_value=0.05; 
% noise amplitude of the traction force data 
% (same unit as traction force data)
% set to zero when unknown
fplot=1; 
% Graphics
% 0: no figures
% 1: plot figures
mult=10^(-2); 
% multiplicative coefficient for a better graphical representation of the
% stress tensor

%% Matrices A and B
 
%% Building 2N*Ninf matrix A such that A*\sigma=T (discretized divergence operator)
% Building matrix Ax (derivative wrt x)
Ax1=diag(ones(C-1, 1), 1)-diag(ones(C, 1), 0); %bulk (order 2)
Ax1=[Ax1,[zeros(C-1, 1); 1]];
Ax=sparse(Ax1);
for i=1:R-1
    Ax=blkdiag(Ax,Ax1); %iteration R times for all the rows
end
clear Ax1;
 
% Building matrix Ay (derivative wrt y)
Ay1=diag(ones((R-1)*C,1), C)-diag(ones(R*C,1 ), 0); %bulk (order 2)
Ay=sparse([Ay1, [zeros((R-1)*C,C); eye(C)]]);
clear Ay1;
 
% Building A from Ax and Ay
A=[Ax/lx, sparse(N,N+C) Ay/ly, sparse(N,N+R); sparse(N,N+R) Ay/ly, sparse(N,N+C), Ax/lx];
clear Ax Ay; 
 
%% Building Ninf*Ninf B covariance matrix of the prior, S_0^{-1}=B/s_0^2
%% Stress norm regularization
B=speye(Ninf);
 
%% Add to the prior a term enforcing the equality of shear components
d1=diag(ones((R-1)*C, 1), C)+diag(ones(R*C,1), 0); %bulk (order 2)
bd1=[d1,[zeros((R-1)*C, C); eye(C)]];
clear d1;
d2=-diag(ones(C-1, 1), 1)-diag(ones(C, 1),0); %bulk (order 2)
d2=[d2,[zeros(C-1, 1); -1]];
bd2=d2;
% loop on rows
for i=1:R-1
   bd2=blkdiag(bd2, d2); 
end
clear d2;
Bdiff=sparse([zeros(N, 2*N+C+R), bd1, bd2]);
clear bd1 bd2;
alpha_xy=10^3; 
% hyperparameter for the equality of shear components (has to be large)
B=B+alpha_xy^2*sparse(Bdiff'*Bdiff);
clear Bdiff;
 
%% Add to the prior a term enforcing the free stress boundary conditions
if BC==1
    Bbcx=zeros(2*C, C*(R+1));
    Bbcy=zeros(2*R, (C+1)*R);   
    Bbcx(1:C,1:C)=speye(C);
    Bbcx(C+1:2*C, 1+R*C:(R+1)*C)=speye(C);
    for i=1:R
         Bbcy(i, 1+(i-1)*(C+1))=1;
         Bbcy(i+R, C+1+(i-1)*(C+1))=1;
    end
    Bbc=sparse([Bbcy, sparse(2*R,(C+1)*R+2*C*(R+1)); sparse(2*C,(C+1)*R), Bbcx, sparse(2*C,(C+1)*R+C*(R+1)); sparse(2*C,(C+1)*R+C*(R+1)), Bbcx, sparse(2*C,(C+1)*R); sparse(2*R,(C+1)*R+2*C*(R+1)), Bbcy]); %xx, yy, yx and xy components
    clear Bbcx Bbcy; 
    % hyperparameter for the free stress boundary condition (has to be large)
    alpha_BC=10^3; 
    B=B+alpha_BC^2*sparse(Bbc'*Bbc);
    clear Bbc;
end
 
%% Traction force
% loop on the time step
for k0=1:kM 

    k0
    
    load(sprintf('%s',ForceName)); % load the traction force data
    f=sprintf('frame%d', k0); % load the k0 time step
    T=traction.(f);
    mTxtrue=T.tx(1:Rtrue,1:Ctrue); % Tx in a R*C matrix form
    mTytrue=T.ty(1:Rtrue,1:Ctrue); % Ty in a R*C matrix form 
    clear T; 
    clear Traction_field;
    [xg, yg]=meshgrid(x,y);
    mTtrue = sqrt(mTxtrue.^2 + mTytrue.^2);

    %% Figure of the traction force field
    if fplot==1
        figure(1000+k0) %field vector T
        quiver(xg,yg,mTxtrue,mTytrue,'b','LineWidth',2);
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('x (\mum)', 'Fontsize', 18)
        ylabel('y (\mum)', 'Fontsize', 18)
        title('t', 'Fontsize', 18)
        axis([xmin-0.01*(xmax-xmin) xmax+0.01*(xmax-xmin) ymin-0.01*(ymax-ymin) ymax+0.01*(ymax-ymin)])
        saveas(1000+k0, ['fig_tractionforcemap_vector_control_'  num2str(nassay, '%d') '_t_' num2str((k0), '%d') '.jpg']  );
        
       figure(5000+k0)
        imagesc(x,y,mTtrue); %ux on rectangular grid
        axis xy
        shading interp
        colorbar
        % caxis([-10 20]);
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('x (\mum)', 'Fontsize', 18)
        ylabel('y (\mum)', 'Fontsize', 18)
        title('Traction force amplitude', 'Fontsize', 18)
        axis([xmin-0.01*(xmax-xmin) xmax+0.01*(xmax-xmin) ymin-0.01*(ymax-ymin) ymax+0.01*(ymax-ymin)])
        saveas(5000+k0, ['fig_tractionforcemap_amplitude_control_'  num2str(nassay, '%d') '_t_' num2str((k0), '%d') '.jpg']  );        
    end
    clear xg yg;
    clear mTtrue;
    
    % coarse-graining
    if nR2 ~= 0 && nC2 ~= 0
        for iC = 1:C
            for iR = 1:R
                mTx(iR,iC) = mean(mean(mTxtrue((iR-1)*nR2+1:(iR+1)*nR2, (iC-1)*nC2+1:(iC+1)*nC2)));
                mTy(iR,iC) = mean(mean(mTytrue((iR-1)*nR2+1:(iR+1)*nR2, (iC-1)*nC2+1:(iC+1)*nC2)));
            end
        end
    elseif nR2 == 0 && nC2 ~= 0
        for iC = 1:C
            for iR = 1:R
                mTx(iR,iC) = mean(mTxtrue(iR, (iC-1)*nC2+1:(iC+1)*nC2));
                mTy(iR,iC) = mean(mTytrue(iR, (iC-1)*nC2+1:(iC+1)*nC2));
            end
        end        
    elseif nR2 ~= 0 && nC2 == 0
        for iC = 1:C
            for iR = 1:R
                mTx(iR,iC) = mean(mTxtrue((iR-1)*nR2+1:(iR+1)*nR2, iC));
                mTy(iR,iC) = mean(mTytrue((iR-1)*nR2+1:(iR+1)*nR2, iC));
            end
        end        
    else 
        mTx = mTxtrue;
        mTy = mTytrue;        
    end
    
    vTx=reshape(mTx', N, 1); % Tx in a N vector form
    vTy=reshape(mTy', N, 1); % Ty in a N vector form 
    T=[vTx; vTy]; % vector T (2N components)
%     [xg, yg]=meshgrid(x,y);
%     
% %% Figure of the traction force field
%     if fplot==1
%         figure(1000+k0) %field vector T
%         quiver(xg,yg,mTx,mTy,'b','LineWidth',2);
%         set(gca, 'FontSize', 18, 'fontName','Times');
%         set(gcf,'Color','w')
%         xlabel('x (\mum)', 'Fontsize', 18)
%         ylabel('y (\mum)', 'Fontsize', 18)
%         title('t', 'Fontsize', 18)
%         axis([xmin-0.01*(xmax-xmin) xmax+0.01*(xmax-xmin) ymin-0.01*(ymax-ymin) ymax+0.01*(ymax-ymin)])
%     end
%     clear xg yg;

%% Hyperparameters

% Lambda
    if meth_Lambda==1
    % MAP estimation with Jeffreys' hyperprior
    % iterative methods 
        % number of steps in the iteration
        stepM=10; 
 
        % initial condition
        step=1;
        Lambda_MAP(step)=10^(-3);
 
        % iteration loop
        while step<stepM
            sigmaMAP=(Lambda_MAP(step)*B+l^2*A'*A)\(l^2*A'*T);
            sigmaMAP_norm(step)=sigmaMAP'*B*sigmaMAP;
            res_norm(step)=(T-A*sigmaMAP)'*(T-A*sigmaMAP);
            step=step+1;
            s0_2=sigmaMAP_norm(step-1)/(4*N+2*(C+R)+2);
            s_2=res_norm(step-1)/(2*N+2);
            Lambda_MAP(step)=l^2*s_2/s0_2;
        end
 
        % graphics
        figure(10+k0)
        semilogy(0:stepM-1,Lambda_MAP,'LineWidth',2)
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('Step','Fontsize',18)
        ylabel('\Lambda', 'Fontsize',18)
        title('MAP estimate', 'Fontsize',18)
 
        % hyperparameter value   
        Lambda=Lambda_MAP(stepM);
        store_Lambda(k0) = Lambda;
        
        % MAP estimate of the noise amplitude
        noise_value_MAP(k0)=sqrt(s_2); 

    elseif meth_Lambda==2
    % (parametric) L-curve
        k=1;
        for lc=-10:0.1:0
             Lambda_Lcurve(k)=10^lc;
             Lsigma_post=(Lambda_Lcurve(k)*B+l^2*A'*A)\(l^2*A'*T);
             sigma_norm(k)=Lsigma_post'*B*Lsigma_post;
             res_norm(k)=(T-A*Lsigma_post)'*(T-A*Lsigma_post);
             k=k+1;
        end
 
        % local slope of the L-curve in logarithmic scale
        dL=(log(sigma_norm(3:end))-log(sigma_norm(1:end-2)))./(log(res_norm(3:end))-log(res_norm(1:end-2)));
        % local curvature of the L-curve
        ddL=2*(dL(3:end)-dL(1:end-2))./(log(res_norm(4:end-1))-log(res_norm(2:end-3)));

        % graphics
        % plot of the L-curve
        figure(100+k0) 
        loglog(res_norm.^0.5,sigma_norm.^0.5,'LineWidth',2)
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('$\mathrm{Residual\,norm\,\Vert\vec{T}-\textbf{A}\vec{\sigma}\Vert\,(kPa)}$','interpreter','latex','Fontsize',18)
        ylabel('Prior norm(kPa.\mum)', 'Fontsize',18)
        title('L curve', 'Fontsize',18)
 
        % plot of the local curvature as a function of Lambda
        figure(200+k0)
        semilogx(Lambda_Lcurve(3:end-2),ddL,'LineWidth',2)
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('$\Lambda$','interpreter','latex','Fontsize',18)
        ylabel('$\mathrm{\frac{d^2\,log(\vec{\sigma}^T\textbf{B}\vec{\sigma})^{1/2}}{d\,log^2\Vert\vec{T}-\textbf{A}\vec{\sigma}\Vert} (\mu m.kPa^{-1})}$', 'interpreter','latex','Fontsize',18) 
        title('Curvature of the L curve', 'Fontsize',18)
     
        % hyperparameter
        % Lambda at the point of maximal curvature
        [max_curv,index_max]=max(ddL);
        Lambda=Lambda_Lcurve(index_max+2); 
        store_Lambda(k0) = Lambda;
        
        % L-curve estimate of the noise amplitude
        noise_value_Lcurve(k0)=sqrt(res_norm(index_max+2)/(2*N+2)); 
        
    elseif meth_Lambda==3
    % set the value of Lambda
        Lambda=1.e-5;
    end


    
%% Stress and covariance matrix estimations
    sigma_post_obs{k0}=(Lambda*B+l^2*A'*A)\(l^2*A'*T);

     if noise_value ~= 0
      Spost{k0}=noise_value^2*l^2*((Lambda*B+l^2*A'*A)\speye(Ninf));
     end
    
end
%clear B;
 
%% Data interpolation and stress representation
for k0=1:kM

    if noise_value ~= 0
        % Inferred variances of xx, yy, xy and yx stress components
        Spost_xx=reshape(diag(Spost{k0}(1:N+R, 1:N+R)), C+1, R)'; 
        Spost_yy=reshape(diag(Spost{k0}(N+R+1:2*N+C+R, N+R+1:2*N+C+R)), C, R+1)'; 
        Spost_xy=reshape(diag(Spost{k0}(2*N+C+R+1:3*N+2*C+R, 2*N+C+R+1:3*N+2*C+R)), C, R+1)'; 
        Spost_yx=reshape(diag(Spost{k0}(3*N+2*C+R+1:4*N+2*(C+R), 3*N+2*C+R+1:4*N+2*(C+R))), C+1, R)'; 
    end        
    vsigma_post_xx=sigma_post_obs{k0}(1:N+R); % inferred stress xx in N+R-vector form
    vsigma_post_yy=sigma_post_obs{k0}(N+R+1:2*N+C+R); % inferred stress yy in N+C-vector form
    vsigma_post_xy=sigma_post_obs{k0}(2*N+C+R+1:3*N+2*C+R); % inferred stress xy in N+C-vector form
    vsigma_post_yx=sigma_post_obs{k0}(3*N+2*C+R+1:4*N+2*(C+R)); % inferred stress yx in N+R-vector form
 
% Interpolation on the grid of the force data
    for i=1:R
        for j=1:C       
            vsigma_xx((i-1)*C+j, 1)=(vsigma_post_xx((i-1)*(C+1)+j)+vsigma_post_xx((i-1)*(C+1)+j+1))/2;
            vsigma_yx((i-1)*C+j, 1)=(vsigma_post_yx((i-1)*(C+1)+j)+vsigma_post_yx((i-1)*(C+1)+j+1))/2;
            vsigma_yy((i-1)*C+j, 1)=(vsigma_post_yy((i-1)*C+j)+vsigma_post_yy((i-1)*C+j+C))/2;
            vsigma_xy((i-1)*C+j, 1)=(vsigma_post_xy((i-1)*C+j)+vsigma_post_xy((i-1)*C+j+C))/2; 
            if noise_value ~= 0
                Ssigma_xx(i,j)=(Spost_xx(i, j)+Spost_xx(i, j+1))/4;
                Ssigma_yx(i,j)=(Spost_yx(i, j)+Spost_yx(i, j+1))/4;
                Ssigma_yy(i,j)=(Spost_yy(i, j)+Spost_yy(i+1, j))/4;
                Ssigma_xy(i,j)=(Spost_xy(i, j)+Spost_xy(i+1, j))/4;
            end
        end
    end
    vsigma_post_xx=vsigma_xx; % inferred stress xx in N-vector form
    vsigma_post_yy=vsigma_yy; % inferred stress yy in N-vector form
    vsigma_post_xy=(vsigma_yx+vsigma_xy)/2; % inferred stress xy in N-vector form
    % Reshape the stress vector in a RxC-matrix form
    sigma_post_xx=reshape(vsigma_post_xx, C, R)'; 
    sigma_post_yy=reshape(vsigma_post_yy, C, R)'; 
    sigma_post_xy=reshape(vsigma_post_xy, C, R)'; 

    if noise_value ~= 0 % estimates of error bars
        Spost_xx=Ssigma_xx;
        Spost_yy=Ssigma_yy;
        Spost_xy=(Ssigma_xy+Ssigma_yx)/4;
        clear Ssigma_xx Ssigma_yy Ssigma_xy Ssigma_yx; % Spost;
        spostxx=Spost_xx.^0.5; 
        spostyy=Spost_yy.^0.5;
        spostxy=Spost_xy.^0.5;
    end
    
    % mean values far from the edges
    scenter(k0) = mean(mean(sigma_post_xy(Rmin:Rmax,Cmin:Cmax)));
%     pcenter(k0) = -( mean(mean(sigma_post_xx(Rmin:Rmax,Cmin:Cmax))) +...
%         mean(mean(sigma_post_xx(Rmin:Rmax,Cmin:Cmax))) )/2.;
    tcenter(k0) = ( mean(mean(sigma_post_xx(Rmin:Rmax,Cmin:Cmax))) +...
        mean(mean(sigma_post_xx(Rmin:Rmax,Cmin:Cmax))) )/2.;
    
%% Figure of the inferred stress tensor field
    if fplot==1
%         for i=1:R
%             for j=1:C
        for i=Rmin:Rmax
            for j=Cmin:Cmax
                % calculation of stress eigenvalues and eigenvectors 
                [V,D]=eig([sigma_post_xx(i,j) sigma_post_xy(i,j); sigma_post_xy(i,j) sigma_post_yy(i,j)]); 
                ab1=[x(j)-mult*D(1,1)*V(1,1),x(j)+mult*D(1,1)*V(1,1)];
                or1=[y(i)-mult*D(1,1)*V(2,1),y(i)+mult*D(1,1)*V(2,1)];
                if D(1,1)>=0
                    % tensile stress in blue
                    figure(2000+k0)
                    hold on
                    plot(ab1,or1,'b','LineWidth',2);
                else
                    % compressive stress in red
                    figure(2000+k0)
                    hold on
                    plot(ab1,or1,'r','LineWidth',2);
                end
                ab2=[x(j)-mult*D(2,2)*V(2,1),x(j)+mult*D(2,2)*V(2,1)];
                or2=[y(i)+mult*D(2,2)*V(1,1),y(i)-mult*D(2,2)*V(1,1)];
                if D(2,2)>=0
                    % tensile stress in blue
                    figure(2000+k0)
                    hold on
                    plot(ab2,or2,'b','LineWidth',2);
                else
                    % compressive stress in red
                    figure(2000+k0)
                    hold on
                    plot(ab2,or2,'r','LineWidth',2);
                end 
            end
        end
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('x (\mum)', 'Fontsize', 18)
        ylabel('y (\mum)', 'Fontsize', 18)
        title('\sigma^{inf}', 'Fontsize', 18)
        axis([xmincenter-0.01*(xmaxcenter-xmincenter) xmaxcenter+0.01*(xmaxcenter-xmincenter) ymincenter-0.01*(ymaxcenter-ymincenter) ymaxcenter+0.01*(ymaxcenter-ymincenter)])
        % axis([xmin-0.01*(xmax-xmin) xmax+0.01*(xmax-xmin) ymin-0.01*(ymax-ymin) ymax+0.01*(ymax-ymin)])
        saveas(2000+k0, ['fig_stresstensor_center_control_' num2str(nassay, '%d') '_t_' num2str((k0), '%d') '.jpg']  );
        
        figure(3000+k0)
        tensioncenter = (sigma_post_xx(Rmin:Rmax,Cmin:Cmax)+sigma_post_yy(Rmin:Rmax,Cmin:Cmax))/2.;
        imagesc(xcenter,ycenter,tensioncenter); %ux on rectangular grid
        axis xy
        shading interp
        colorbar
        % caxis([-10 20]);
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('x (\mum)', 'Fontsize', 18)
        ylabel('y (\mum)', 'Fontsize', 18)
        title('tension', 'Fontsize', 18)
        axis([xmincenter-0.01*(xmaxcenter-xmincenter) xmaxcenter+0.01*(xmaxcenter-xmincenter) ymincenter-0.01*(ymaxcenter-ymincenter) ymaxcenter+0.01*(ymaxcenter-ymincenter)])
        % axis([xmin-0.01*(xmax-xmin) xmax+0.01*(xmax-xmin) ymin-0.01*(ymax-ymin) ymax+0.01*(ymax-ymin)])
        saveas(3000+k0, ['fig_tensionmap_center_control_'  num2str(nassay, '%d') '_t_' num2str((k0), '%d') '.jpg']  );
    end
 
%% Accuracy of the inference
    Tinf=A*sigma_post_obs{k0}; 
    % divergence of the inferred stress Tinf =A*sigma_inf
    clear sigma_post_obs{k0};
    Tinfx=Tinf(1:N); 
    Tinfy=Tinf(N+1:2*N);
    clear Tinf;
    R2x=1-((vTx'-Tinfx')*(vTx-Tinfx))/((vTx'-mean(vTx))*(vTx-mean(vTx)));
    R2y=1-((vTy'-Tinfy')*(vTy-Tinfy))/((vTy'-mean(vTy))*(vTy-mean(vTy)));
    % R2: coefficient of determination between the data and the traction forces 
    % obtained by calculating the divergence of the inferred stress 
    % (should be close to 1)
    R2=(R2x+R2y)/2; 
    store_R2(k0) = R2;

% Inferred stress mean values
    sxxm=sum(vsigma_post_xx)/N;
    syym=sum(vsigma_post_yy)/N;
    sxym=sum(vsigma_post_xy)/N;
 
% Inferred hyperparameter
    stress.Lambda{k0}=Lambda;
    
% Inferred stress values
    stress.sxx{k0}=sigma_post_xx;
    stress.syy{k0}=sigma_post_yy;
    stress.sxy{k0}=sigma_post_xy;
    stress.x{k0}=x;
    stress.y{k0}=y;
    stress.R2_T{k0}=R2;
    stress.mean_from_sigma{k0}=[sxxm,syym,sxym];   

% Stress mean values from traction force
    if BC==1
        sxx_mean=-sum(sum(reshape(vTx, C, R)'.*(ones(R, 1)*(x-(xmax-xmin)/2))))/N;
        syy_mean=-sum(sum(reshape(vTy, C, R)'.*((y-(ymax-ymin)/2)'*ones(1, C))))/N;
        sxy_mean=-sum(sum(reshape(vTy, C, R)'.*(ones(R, 1)*(x-(xmax-xmin)/2))+reshape(vTx, C, R)'.*((y-(ymax-ymin)/2)'*ones(1,C))))/(2*N);
        stress.mean_from_t{k0}=[sxx_mean,syy_mean,sxy_mean];
    end
    
% Values of error bars on the inferred stress
    if noise_value~=0
        stress.error_sxx{k0}=spostxx;
        stress.error_syy{k0}=spostyy;
        stress.error_sxy{k0}=spostxy;
    end
end

% save output
% save(['Stress_BISM_Before_Caly_' num2str((k0), '%d') '.mat'] ,'stress');
save(['Stress_BISM_control_' num2str(nassay, '%d') '_t_' num2str((k0), '%d') '.mat'] ,'stress');



figure(1013)
plot(1:kM,tcenter,'+');
set(gca, 'FontSize', 18, 'fontName','Times');
set(gcf,'Color','w')
xlabel('frame', 'Fontsize', 18)
ylabel('central tension', 'Fontsize', 18)
% saveas(1013, ['fig_tractionforce_pcenter_Before_Caly_'  num2str((k0), '%d') '.jpg']  );
saveas(1013, ['fig_tractionforce_tcenter_control_' num2str(nassay, '%d') '_t_' num2str((k0), '%d') '.jpg']  );


figure(1012)
plot(1:kM,scenter,'+');
set(gca, 'FontSize', 18, 'fontName','Times');
set(gcf,'Color','w')
xlabel('frame', 'Fontsize', 18)
ylabel('central shear stress', 'Fontsize', 18)
% saveas(1012, ['fig_tractionforce_scenter_Before_Caly_' num2str((k0), '%d') '.jpg']  );
saveas(1012, ['fig_tractionforce_scenter_control_' num2str(nassay, '%d') '_t_' num2str((k0), '%d') '.jpg']  );



        
        
