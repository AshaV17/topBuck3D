%The code is developed from original code topBuck250.m developed in the paper - F. Ferrari, O. Sigmund, and J. K. Guest, “Topology optimization with linearized buckling criteria in 250 lines of Matlab,” Struct Multidisc Optim, vol. 63, no. 6, pp. 3045–3066, Jun. 2021, doi: 10.1007/s00158-021-02854-x.
%main code
function topBuck3D(nelx,nely,nelz,penalK,rmin,ft,ftBC,eta,beta,maxit,ocPar,Lx,Ly,penalG,nEig,pAgg,prSel,fil)
% ---------------------------- PRE. 1) MATERIAL AND CONTINUATION PARAMETERS
[E0,Emin,nu] = deal(1,1e-9,0.3);                                           % Young's moduli & Poisson's ratio
penalCntK = {25,1,25,0.25};                                                 % continuation scheme on K-penal
penalCntG = {25,1,25,0.25};                                                 % " " on G-penal
betaCnt  = {150,18,25,2};                                                     % " " on beta
pAggCnt  = {2e5,1,25,2};                                                   % " " on the KS aggregation factor
if prSel{1}(1) == 'V', volfrac = 1.0; else, volfrac = prSel{2}(end); end   % initialize volume fraction
% ----------------------------------------- PRE. 2) DISCRETIZATION FEATURES
Lz = nelz/nelx*Lx;                                                         % recover Lz from aspect ratio#3D#
nEl = nelx * nely * nelz;                                                  % number of elements          #3D#
elNrs = reshape(1:nEl,nely,nelz,nelx);                                          % element numbering
nodeNrs = int32( reshape( 1 : ( 1 + nelx ) * ( 1 + nely ) * ( 1 + nelz ), ...
    1 + nely, 1 + nelz, 1 + nelx ) );                                      % nodes numbering             #3D#
cVec = reshape( 3 * nodeNrs( 1 : nely, 1 : nelz, 1 : nelx ) + 1, nEl, 1 ); %                             #3D#
cMat = cVec+int32( [0,1,2,3*(nely+1)*(nelz+1)+[0,1,2,-3,-2,-1],-3,-2,-1,3*...
   (nely+1)+[0,1,2],3*(nely+1)*(nelz+2)+[0,1,2,-3,-2,-1],3*(nely+1)+...
   [-3,-2,-1]]);                                                           % connectivity matrix         #3D#
nDof = ( 1 + nely ) * ( 1 + nelz ) * ( 1 + nelx ) * 3;                     % total number of DOFs        #3D#
if ftBC == 'N', bcF = 'symmetric'; else, bcF = 0; end                      % filter BC selector
% ---------------------------------------------- elemental stiffness matrix
[ sI, sII ] = deal( [ ] );
for j = 1 : 24
    sI = cat( 2, sI, j : 24 );
    sII = cat( 2, sII, repmat( j, 1, 24 - j + 1 ) );
end
[ iK , jK ] = deal( cMat( :,  sI )', cMat( :, sII )' );
Iar = sort( [ iK( : ), jK( : ) ], 2, 'descend' );                          % reduced assembly indexing
Ke = 1/(1+nu)/(2*nu-1)/144 *( [ -32;-6;-6;8;6;6;10;6;3;-4;-6;-3;-4;-3;-6;10;...
    3;6;8;3;3;4;-3;-3; -32;-6;-6;-4;-3;6;10;3;6;8;6;-3;-4;-6;-3;4;-3;3;8;3;...
    3;10;6;-32;-6;-3;-4;-3;-3;4;-3;-6;-4;6;6;8;6;3;10;3;3;8;3;6;10;-32;6;6;...
    -4;6;3;10;-6;-3;10;-3;-6;-4;3;6;4;3;3;8;-3;-3;-32;-6;-6;8;6;-6;10;3;3;4;...
    -3;3;-4;-6;-3;10;6;-3;8;3;-32;3;-6;-4;3;-3;4;-6;3;10;-6;6;8;-3;6;10;-3;...
    3;8;-32;-6;6;8;6;-6;8;3;-3;4;-3;3;-4;-3;6;10;3;-6;-32;6;-6;-4;3;3;8;-3;...
    3;10;-6;-3;-4;6;-3;4;3;-32;6;3;-4;-3;-3;8;-3;-6;10;-6;-6;8;-6;-3;10;-32;...
    6;-6;4;3;-3;8;-3;3;10;-3;6;-4;3;-6;-32;6;-3;10;-6;-3;8;-3;3;4;3;3;-4;6;...
    -32;3;-6;10;3;-3;8;6;-3;10;6;-6;8;-32;-6;6;8;6;-6;10;6;-3;-4;-6;3;-32;6;...
    -6;-4;3;6;10;-3;6;8;-6;-32;6;3;-4;3;3;4;3;6;-4;-32;6;-6;-4;6;-3;10;-6;3;...
    -32;6;-6;8;-6;-6;10;-3;-32;-3;6;-4;-3;3;4;-32;-6;-6;8;6;6;-32;-6;-6;-4;...
    -3;-32;-6;-3;-4;-32;6;6;-32;-6;-32]+nu*[ 48;0;0;0;-24;-24;-12;0;-12;0;...
    24;0;0;0;24;-12;-12;0;-12;0;0;-12;12;12;48;0;24;0;0;0;-12;-12;-24;0;-24;...
    0;0;24;12;-12;12;0;-12;0;-12;-12;0;48;24;0;0;12;12;-12;0;24;0;-24;-24;0;...
    0;-12;-12;0;0;-12;-12;0;-12;48;0;0;0;-24;0;-12;0;12;-12;12;0;0;0;-24;...
    -12;-12;-12;-12;0;0;48;0;24;0;-24;0;-12;-12;-12;-12;12;0;0;24;12;-12;0;...
    0;-12;0;48;0;24;0;-12;12;-12;0;-12;-12;24;-24;0;12;0;-12;0;0;-12;48;0;0;...
    0;-24;24;-12;0;0;-12;12;-12;0;0;-24;-12;-12;0;48;0;24;0;0;0;-12;0;-12;...
    -12;0;0;0;-24;12;-12;-12;48;-24;0;0;0;0;-12;12;0;-12;24;24;0;0;12;-12;...
    48;0;0;-12;-12;12;-12;0;0;-12;12;0;0;0;24;48;0;12;-12;0;0;-12;0;-12;-12;...
    -12;0;0;-24;48;-12;0;-12;0;0;-12;0;12;-12;-24;24;0;48;0;0;0;-24;24;-12;...
    0;12;0;24;0;48;0;24;0;0;0;-12;12;-24;0;24;48;-24;0;0;-12;-12;-12;0;-24;...
    0;48;0;0;0;-24;0;-12;0;-12;48;0;24;0;24;0;-12;12;48;0;-24;0;12;-12;-12;...
    48;0;0;0;-24;-24;48;0;24;0;0;48;24;0;0;48;0;0;48;0;48 ] );             % elemental stiffness matrix  #3D#
Ke0( tril( ones( 24 ) ) == 1 ) = Ke';
Ke0 = reshape( Ke0, 24, 24 );
Ke0 = Ke0 + Ke0' - diag( diag( Ke0 ) );                                    % recover full matrix
if any(prSel{1}=='B') % >>>>>>>>>>>>>>>> PERFORM ONLY IF BUCKLING IS ACTIVE #B#
    Cmat0 = [1-nu,nu,nu,0,0,0;nu,1-nu,nu,0,0,0;nu,nu,1-nu,0,0,0;...
        0,0,0,(1-2*nu)/2,0,0;0,0,0,0,(1-2*nu)/2,0;0,0,0,0,0,(1-2*nu)/2]...
        /(1+nu)/(1-2*nu);                                                  % non-dimensional elasticity matrix  #3D#
    xiG = sqrt(1/3)*[-1,1]; etaG = xiG;psiG = xiG; wxi = [1,1]; weta = wxi;% Gauss nodes and weights #3D#
    wpsi = wxi;
    xe = [-1,-1,-1;1,-1,-1;1,1,-1;-1,1,-1;-1,-1,1;1,-1,1;1,1,1;-1,1,1].*Lx/nelx; % dimensions of the elements #3D#
    lMat = zeros(6, 9); lMat(1, 1) = 1; lMat(2, 5) = 1; lMat(3, 9) = 1;    % placement matrix #3D#
    lMat(4, 2) = 1; lMat(4, 4) = 1; lMat(5, 3) = 1; lMat(5, 7) = 1;  
    lMat(6, 6) = 1; lMat(6, 8) = 1;
    dN = @(xi,zi,yi) 0.125*[(zi-1)*(1-yi),(1-zi)*(1-yi),(1+zi)*(1-yi),...
        (-1-zi)*(1-yi),(zi-1)*(1+yi),(1-zi)*(1+yi),(1+zi)*(1+yi),(-1-zi)*...
        (1+yi); (xi-1)*(1-yi),(-1-xi)*(1-yi),(1+xi)*(1-yi),(1-xi)*(1-yi),...
        (xi-1)*(1+yi),(-1-xi)*(1+yi),(1+xi)*(1+yi),(1-xi)*(1+yi);...
        (1-xi)*(zi-1),(xi+1)*(zi-1),(-1-xi)*(1+zi),(1-xi)*(-1-zi),...
        (1-xi)*(1-zi),(1+xi)*(1-zi),(1+xi)*(1+zi),(1-xi)*(1+zi)];          % shape funct. logical derivatives #3D#
    B0 = @(gradN) lMat * kron(gradN,eye(3));                               % strain-displacement matrix
    [indM,t2ind] = deal([1,4,7,10,13,16,19,22,70,73,76,79,82,85,88,130,133,...
        136,139,142,145,181,184,187,190,193,223,226,229,232,256,259,262,...
        280,283,295],[2,3,4,5,6,7,8,10,11,12,13,14,15,17,18,19,20,21,23,24, ...
        25,26,28,29,30,32,33,35]);                                         % auxiliary set of indices (1) #3D#
    [iG,jG] = deal(iK(indM,:),jK(indM,:));                                 % indexing of unique G coefficients
    IkG = sort([iG(:), jG(:)],2,'descend'); clear iG jG iK jK;                   % indexing G entries (lower half)
    [a1,a2]=deal(reshape(IkG(:,2),36,nEl)', reshape(IkG(:,1),36,nEl)');    % auxiliary set of indices (2) #3D#
    dZdu = zeros(36,24);                                                   % build U-derivative of matrix G #3D#
    for ii = 1 : 24                                                        % loop on the displacement components #3D#
        tt = 0; Uvec = zeros(24,1); Uvec(ii,1) = 1;                        % set a single displ. component #3D#
        se = Cmat0*B0((dN(0,0,0)*xe)\dN(0,0,0))*Uvec;                      % stresses at the element center #3D#
        for j = 1 : length(xiG)
            for k = 1 : length(etaG)
                for kk = 1 : length(psiG)
                    xi = xiG(j); zi = etaG(k); yi= psiG(kk);               % current integration points #3D#
                    w = wxi(j)*weta(k)*wpsi(kk)*det(dN(xi,zi,yi)*xe);      % current integration weight #3D#
                    gradN = (dN(xi,zi,yi)*xe)\dN(xi,zi,yi);                % shape funct. physical derivatives #3D#
                    B1 = [kron(gradN,[1,0,0]); kron(gradN,[0,1,0]);...
                        kron(gradN,[0,0,1])];                              % deformation gradient #3D#
                    tt = tt+(B1'*kron(eye(3),[se(1),se(4),se(5);...
                        se(4),se(2),se(6);se(5),se(6),se(3)])*B1)*w;                                                % current contribution to dG/du_i #3D#
                end
            end
        end
        dZdu(:,ii) = tt([1,4,7,10,13,16,19,22,76,79,82,85,88,91,94,151,154,...
        157,160,163,166,226,229,232,235,238,301,304,307,310,376,379,382,...
        451,454,526])';                                                    % extract independent coefficients #3D#
    end
    dZdu(t2ind,:) = 2*dZdu(t2ind,:);                                       % x2 columns for v-m-v product
    fKS = @(p,v)max(v)+log(sum(exp(p*(v-max(v)))))/p;                      % KS aggregation function
    dKS = @(p,v,dv)sum(exp(p.*(v-max(v)))'.*dv,2)./sum(exp(p.*(v-max(v))));% derivative of the KS function
end % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #B#
% ----------------------------- PRE. 3) LOADS, SUPPORTS AND PASSIVE DOMAINS
% to be changed for different structures and cell sizes. This is for SC lattice of 20 x 20 x 20 size
tkx=4;
tky=4;
tkz=4;
tkxb=nelx/10;
tkyb=nely/10;
tkzb=nelz/10;
lcDof = unique(cMat(unique([elNrs(: , [1:tkz,end-tkz+1:end],end), ...
    (elNrs([1:tky,end-tky+1:end] , :,end))']),[4,7,16,19]));%node numbers on loading is 2,3,6,7 x direction dofs only
lcDofix = unique(cMat(unique([elNrs(: , [1:tkz,end-tkz+1:end],end), ...
    (elNrs([1:tky,end-tky+1:end] , :,end))']),[5,6,8,9,17,18,20,21]));%node numbers on loading is 2,3,6,7 y,z direction dofs only
fixed = unique(cMat(unique([elNrs(: ,[1:tkz,end-tkz+1:end],1),...     %node numbers on BC is 1,4,5,8 all dofs
     (elNrs([1:tky,end-tky+1:end] , :,1))']),[1:3,10:15,22:24]));
edgedof=unique(cMat(unique([elNrs(: , [1,end],end); (elNrs([1,end], : ,end))']),[4,7,16,19]));
modF = 2e-3;                                                                % modulus of the force density
F = fsparse(lcDof(:),1,-modF,[nDof,1]);                                     % define load vector
F(edgedof(:)) = F(lcDof(1))/2;                                            % consistent load on end nodes

pa1 = elNrs([1,end], [1,end] ,:);                        % outer edge
pa2 = elNrs([1,end], [tkz,end-tkz+1] ,:);                % outer and insides
pa3 = elNrs([tky,end-tky+1], [1,end] ,:);                % outer and insides
pa4 = elNrs([1:2,tky-1:tky,end-tky+1:end-tky+2,end-1:end], [1:2,tkz-1:tkz,end-tkz+1:end-tkz+2,end-1:end] ,:);                          % inner edge beams
pa5 = elNrs(:, [1:tkz,end-tkz+1:end], [1:tkxb,end-tkxb+1:end]);%top bottom beams 
pa6 = elNrs([1:tky,end-tky+1:end],:, [1:tkxb,end-tkxb+1:end]);%front back beams 
%center and faces void
pa7 = elNrs(tky+1:end-tky, tkz+1:end-tkz, [1:tkxb end-tkxb:end tkxb+1:end-tkxb]);                                              % interior void vol
pa8 = elNrs(tky+1:end-tky, [1:tkz end-tkz:end], tkxb+1:end-tkxb);
pa9 = elNrs([1:tky end-tky:end], tkz+1:end-tkz, tkxb+1:end-tkxb);

[pasS,pasV] = deal(unique([pa1(:);pa2(:);pa3(:);pa5(:);pa6(:)]),unique([... 
    pa7(:);pa8(:);pa9(:)]));                     % define passive domains

free = setdiff(1:nDof, unique([fixed;lcDofix]));                                             % set of free DOFs
act = setdiff((1:nEl)',union(pasS(:),pasV(:)));                            % set of active design variables

% --------------------------------------- PRE. 4) DEFINE IMPLICIT FUNCTIONS
prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v(:)-eta)))./...
    (tanh(beta*eta)+tanh(beta*(1-eta)));                                   % projection
deta = @(v,eta,beta) - beta * csch( beta ) .* sech( beta * ( v( : ) - eta ) ).^2 .* ...
    sinh( v( : ) * beta ) .* sinh( ( 1 - v( : ) ) * beta );                % projection eta-derivative 
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));% proj. x-derivative
cnt = @(v,vCnt,l) v+(l>=vCnt{1}).*(v<vCnt{2}).*(mod(l,vCnt{3})==0).*vCnt{4};
% -------------------------------------------------- PRE. 5) PREPARE FILTER
[dy,dz,dx]=meshgrid(-ceil(rmin)+1:ceil(rmin)-1,...
    -ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1 );
h = max( 0, rmin - sqrt( dx.^2 + dy.^2 + dz.^2 ) );                        % conv. kernel                #3D#
Hs = imfilter( ones( nely, nelz, nelx ), h, bcF );                         % matrix of weights (filter)  #3D#
dHs = Hs;
% ------------------------ PRE. 6) ALLOCATE AND INITIALIZE OTHER PARAMETERS
[x,dsK,dsG,dmKS,dV] = deal(zeros(nEl,1));                                  % initialize vectors of size nElx1
U = zeros(nDof,1);                                                           % " " of size nDofx1 
dV(act,1) = 1/nEl/volfrac;                                                         % derivative of volume fraction
[xpOld,loop,restartAs,ch,plotL,plotR,muVec] = deal(0,0,0,1,[],[],[]);      % misc array & parameters
if nargin > 17
  load(fil); x = xInitial;                                   % initialize design from saved data
else
    x(act) = volfrac;%(volfrac*(nEl-length(pasV))-length(pasS))/length(act);        % volume fraction on "active" set
    x(pasS) = 1;                                                           % set x=1 on "passive solid" set
end
xPhys = x; clear iG jG dx dy;                                        % initialize xPhys and free memory
% ================================================= START OPTIMIZATION LOOP
while ch > 1e-6 && loop < maxit
  loop = loop + 1;                                                         % update iter. counter
  % ----------- RL. 1) COMPUTE PHYSICAL DENSITY FIELD (AND ETA IF PROJECT.)
  xTilde = imfilter( reshape( x, nely, nelz, nelx ), h, bcF ) ./ Hs;       % filtered field              #3D#
  xPhys( act ) = xTilde( act );                                            % reshape to column vector
  if ft > 1                              % compute optimal eta* with Newton
      f = ( mean( prj( xPhys, eta, beta ) ) - volfrac )  * (ft == 3);      % function (volume)
      while abs( f ) > 1e-6  && prSel{1}(1) ~= 'V'         % Newton process for finding opt. eta
          eta = eta - f / mean( deta( xPhys, eta, beta ) );
          f = mean( prj( xPhys, eta, beta ) ) - volfrac;
      end
      dHs = Hs ./ reshape( dprj( xPhys, eta, beta ), nely, nelz, nelx );   % sensitivity modification    #3D#
      xPhys = prj( xPhys, eta, beta );                                     % projected (physical) field
  end
  ch = max(abs(xPhys-xpOld));
  xpOld = xPhys;
  % -------------------------- RL. 2) SETUP AND SOLVE EQUILIBRIUM EQUATIONS
  sK = (Emin+xPhys.^penalK*(E0-Emin));                                     % stiffness interpolation
  dsK(act) = penalK*(E0-Emin)*xPhys(act).^(penalK-1);                      % derivative of " "
  sK = reshape(Ke(:)*sK',length(Ke)*nEl,1);
  K = fsparse(Iar(:,1),Iar(:,2),sK,[nDof,nDof]);                           % assemble stiffness matrix
  K = K+K'-diag(diag(K));                                                  % symmetrization of K
  dK = decomposition(K(free,free),'chol','lower');                         % decompose K and store factor
  U(free) = dK \ F(free);                                                  % solve equilibrium system
  dc = -dsK.*sum((U(cMat)*Ke0).*U(cMat),2);                                % compute compliance sensitivity
  % ---------------------- RL. 4) SELECT OBJECTIVE FUNCTION AND CONSTRAINTS
  if loop==1, c0=F'*U; v0=mean(xPhys(:)); end % initial compliance & volume fraction
  switch prSel{1}                    % select optimization problem to solve
      case ['C','V']           % minimize compliance with volume constraint
          g0 = F'*U/c0;
          dg0 = imfilter(reshape(dc/c0,nely,nelz,nelx)./dHs,h,bcF);        %  #3D
          g1 = mean(xPhys(:))/volfrac-1;
          dg1 = imfilter( reshape( dV/volfrac, nely, nelz, nelx ) ./ dHs, h, bcF );   %  #3D
      case ['V','C']           % minimize volume with compliance constraint
          g0 = mean(xPhys(:))./v0;
          dg0 = imfilter(reshape(dV/v0,nely,nelz,nelx)./dHs,h,bcF);        %  #3D
          g1 = (F'*U)/(prSel{2}*c0)-1;
          dg1 = imfilter(reshape(dc/(prSel{2}*c0),nely,nelz,nelx)./dHs,h,bcF);%  #3D
      case ['B','C','V']% maximize BLF with compliance & volume constraints (Eq. 13 paper)
          if loop==1, muKS0=fKS(pAgg,mu(1:nEig)); g0=1; cMax=prSel{2}(1);
          else, g0=fKS(pAgg,mu(1:nEig))/muKS0; end                         % KS aggregation of mu (=1/lambda)
          dmKS = dKS(pAgg,mu(1:nEig),dmu);                                 % KS aggregation of dmu
          dg0 = imfilter(reshape(dmKS/muKS0,nely,nelz,nelx)./dHs,h,bcF);   % back-filter KS sensitivity  %  #3D
          % -- Constraint function: KS aggregation of compliance and volume
          g1Vec = [F'*U;mean(xPhys(:))]./[cMax*c0;volfrac]-1;              % set of constraints ['C','V']
          dg1c = imfilter(reshape(dc/(cMax*c0),nely,nelz,nelx)./dHs,h,bcF);% back-filter compliance derivative  %  #3D
          dg1V =imfilter( reshape( dV/volfrac, nely, nelz, nelx ) ./ dHs, h, bcF );  % back-filter volume derivative  %  #3D
          g1 = fKS(pAgg,g1Vec);                                            % aggregate the two constraints
          dg1 = dKS(pAgg,g1Vec,[dg1c(:),dg1V(:)]);                         % sensitivity of the KS constraint
          plotL(loop,:) = [1/g0/muKS0,1/mu(1)]; strL='KS(-),\lambda_1(--)';
          plotR(loop,:) = [g1,g1Vec']; strR='g_1(-),gC(--),gV(.-)';
          muVec(loop,:) = mu';
      case ['V','C','B']    % min volume with compliance & BLF constraints (Eq. 14 paper)
          g0 = mean(xPhys(:))./v0;
          dg0 = imfilter( reshape( dV/volfrac, nely, nelz, nelx ) ./ dHs, h, bcF );   %  #3D
          % ---- Constraint function: KS aggregation of BLFs and compliance
          muKS = fKS(pAgg,mu(1:nEig));                                     % KS aggregation of mu
          dmKS = dKS(pAgg,mu(1:nEig),dmu);                                 % KS aggregation of dmu
          g1Vec = [prSel{2}(2)*muKS;F'*U]./[1;prSel{2}(1)*c0]-1;           % set of constraints 'B','C'
          dg1l = imfilter(reshape(dmKS*prSel{2}(2),nely,nelz,nelx)./dHs,h,bcF); % back-filter dmu  %  #3D
          dg1c = imfilter(reshape(dc/(prSel{2}(1)*c0),nely,nelz,nelx)./dHs,h,bcF);% back-filter dc  %  #3D
          g1 = fKS(pAgg,g1Vec);                                            % aggregate the two constraints
          dg1 = dKS(pAgg,g1Vec,[dg1l(:),dg1c(:)]);                         % sensitivity of the KS constraint
          plotL(loop,:) = g0; strL = 'g_0';
          plotR(loop,:) = [g1,g1Vec']; strR='g_1(-),gL(--),gC(.-)';
          muVec = cat(1,muVec,mu');
  end

  %----------------- RL. 4) UPDATE DESIGN VARIABLES AND APPLY CONTINUATION
  if loop==1, xOld = x(act); xOld1 = xOld; as = []; end                    % initialize MMA history parameters
  [xo,as,lmid]=ocUpdate(loop,x(act),dg0(act),g1,dg1(act),ocPar,xOld,xOld1,as,beta,restartAs);
  xOld1 = xOld; xOld = x(act); x(act) = xo;
  %  apply continuation on penalization(s), beta & aggregation parameter(s)
  penalKold = penalK; penalGold = penalG; betaOld = beta;
  [penalK,penalG,beta] = deal(cnt(penalK, penalCntK, loop), ...
      cnt(penalG,penalCntG,loop), cnt(beta,betaCnt,loop) );
  if (beta-betaOld~= 0 || penalK-penalKold~=0 || penalG-penalGold~=0)
            restartAs = 1; else, restartAs = 0; end                              % restart asymptotes if needed

  % ----------------------------------------- RL. 8) PRINT AND PLOT RESULTS 
  fprintf('It.:%2i g0:%7.4f g1:%0.2e penalK:%7.2f penalG:%7.2f eta:%7.2f beta:%7.1f ch:%0.3e lm:%0.3e\n', ...
    loop,g0,g1,penalK,penalG,eta,beta,ch,lmid);
  isovals = shiftdim( reshape( xPhys, nely, nelz, nelx ), 2 );              %  #3D#
  isovals = smooth3( isovals, 'box', 1 );
  patch(isosurface(isovals, .5),'FaceColor','b','EdgeColor','none');
  patch(isocaps(isovals, .5),'FaceColor','r','EdgeColor','none');
  drawnow; view( [ 145, 25 ] ); axis equal tight off; %cla();
end
end
%%
function [ x, as, lmid ] = ocUpdate( loop, xT, dg0, g1, dg1, ocPar, xOld, xOld1, ...
    as, beta, restartAs )
% -------------------------------- definition of asymptotes and move limits
[xU,xL] = deal(min(xT+ocPar(1),1), max(xT-ocPar(1),0));
if (loop<2.5 || restartAs==1)
    as = xT+[-0.5,0.5].*(xU-xL)./(beta+1);
else
    tmp = (xT-xOld).*(xOld-xOld1);
    gm = ones(length(xT),1);
    [gm(tmp>0), gm(tmp<0)] = deal(ocPar(3),ocPar(2));
    as = xT+gm.*[-(xOld-as(:,1)),(as(:,2)-xOld)];
end
xL = max( 0.9*as(:,1)+0.1*xT,xL);                    % adaptive lower bound
xU = min( 0.9*as(:,2)+0.1*xT,xU);                    % adaptive upper bound
% ----- split (+) and (-) parts of the objective and constraint derivatives
p0_0 = (dg0>0).*dg0; q0_0 = (dg0<0).*dg0;
p1_0 = (dg1>0).*dg1; q1_0 = (dg1<0).*dg1;
[p0,q0] = deal(p0_0.*(as(:,2)-xT).^2,-q0_0.*(xT-as(:,1)).^2);
[p1,q1] = deal(p1_0.*(as(:,2)-xT).^2,-q1_0.*(xT-as(:,1)).^2);
% ---------------------------------------- define the primal projection map
primalProj = @(lm) min(xU,max(xL,(sqrt(p0+lm*p1).*as(:,1)+sqrt(q0+lm*q1).*as(:,2))...
    ./(sqrt(p0+lm*p1)+sqrt(q0+lm*q1))));
psiDual = @(lm) g1 - ( (as(:,2)-xT)'*p1_0 - (xT-as(:,1))'*q1_0 ) + ...
    sum(p1./(max(as(:,2)-primalProj(lm),1e-12)) + q1./(max(primalProj(lm)-as(:,1),1e-12)));
% ------------- compute the Lagrange multiplier and update design variables
lmUp = 1e6; x = xT; lmid = -1;
if psiDual( 0 ) * psiDual( lmUp ) < 0  % check if LM is within the interval
    lmid = fzero( psiDual, [ 0, lmUp ] );
    x = primalProj( lmid );
elseif psiDual( 0 ) < 0                       % constraint cannot be active
   lmid = 0;
   x = primalProj( lmid );
elseif psiDual( lmUp ) > 0                 % constraint cannot be fulfilled
   lmid = lmUp;
   x = primalProj( lmid );
end
%
end
% to call function, use topBuck3D(20,20,20,3,1.5,2,'N',0.5,2,500,[0.1,0.7,1.2],1,1,3,12,200,{['B','C','V'],[2.5,0.15]},'ig.mat') 
% after running topBuck3D(20,20,20,3,1.5,2,'N',0.5,2,500,[0.1,0.7,1.2],1,1,[],[],[],{['V','C'],2.5}) and saving densities in 'ig.mat' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was originally written by F. Ferrari, O. Sigmund        %
% Dept. of Solid Mechanics-Technical University of Denmark,2800 Lyngby (DK)%
% edited by A. Viswanath of Khalifa University for 3D changes.             %
% Please send your comments to: asha.viswanath@gmail.com                   %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper Ferrari, F. Sigmund, O. - A new generation 99 %
% line Matlab code for compliance Topology Optimization and its extension  %
% to 3D, SMO, 2020    and                                                  %
% Asha Viswanath, Wesley James Cantwell, Kamran Ahmed Khan, Enhancing      %
% Buckling Resistance of Strut Lattice Structures via Three-Dimensional    % 
% Topology Optimization, Journal of Computational Design and Engineering,  % 
% 2025                                                                     %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
