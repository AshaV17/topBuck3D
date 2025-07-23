% 1.	SC Structure (20mm,15% relative density)
% -------------------------- PRE. 3) LOADS, SUPPORTS AND PASSIVE DOMAINS
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
modF = 2e-3;

% column borders
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
free = setdiff(1:nDof, unique([fixed;lcDofix])); % set of free DOFs

% 2.	BCC Structure (40mm,10% relative density)
% -------------------------- PRE. 3) LOADS, SUPPORTS AND PASSIVE DOMAINS
[voxel,Density] = GenerateVoxel(40,'bcc.txt',0.15);
pasV=find(voxel==0);
tkx=5;
tky=5;
tkz=5;
lc =elNrs([1:tky,end-tky+1:end], [1:tkz,end-tkz+1:end],end );
fixe = elNrs([1:tky,end-tky+1:end] , [1:tkz,end-tkz+1:end],1);
lcDof = unique(cMat(lc,[4,7,16,19]));  %node numbers on loading is 2,3,6,7 x direction dofs only
lcDofix = unique(cMat(lc,[5,6,8,9,17,18,20,21])); %node numbers on loading is 2,3,6,7 x direction dofs only
fixed = unique(cMat(fixe,[1:3,10:15,22:24]));
modF = 2e-3;
edgedof=unique(cMat(unique([[elNrs([1,tkz,end-tkz+1,end], ...
    [1:tkz,end-tkz+1:end],end )],[elNrs([1:tkz,end-tkz+1:end],...
    [1,tkz,end-tkz+1,end], end )]']),[4,7,16,19]));
pasS = [lc(:);fixe(:)];

% 3.	nUnit x nUnit x 1 Structure
% ----------------------------- PRE. 3) LOADS, SUPPORTS AND PASSIVE DOMAINS
tkx=4;
tky=4;
tkz=4;
tkxb=2;
tkyb=2;
tkzb=2;
%general indices for n unit lattices
coordX=[];
coordY=[];
ctx=[];
cty=[];
ctx1=[];
cty1=[];
cvy=[];
cvy1=[];
cvxb=[];
cvyb=[];
cvxf=[];
cvyf=[];
for uC=1:nUnit-1
    coordX(uC)=uC*nelxU;
    coordY(uC)=uC*nelyU;
    ctx=[ctx coordX(uC)-tkxb+1:coordX(uC)];
    cty=[cty coordY(uC)-tkyb+1:coordY(uC)];
    ctx1=[ctx1 coordX(uC)+1:coordX(uC)+tkxb];
    cty1=[cty1 coordY(uC)+1:coordY(uC)+tkyb];
    cvy=[cvy (uC-1)*(nelyU-1)+tky+1:uC*nelyU-tky uC*nelyU+tky:(uC+1)*nelyU-tky+1];
    cvyf=[cvyf (uC-1)*(nelyU-1)+tky+1:uC*nelyU-tky uC*nelyU+tky+1:(uC+1)*nelyU-tky];
    cvxf=[cvxf (uC-1)*(nelxU-1)+tkxb+1:nelxU-tkxb uC*nelxU+tkxb:(uC+1)*nelxU-tkxb];
    cvy1=[cvy1 coordY(uC)-tkx:coordY(uC)+tky];
    cvyb=[cvyb (uC-1)*(nelyU-1)+tkyb+1:uC*nelyU uC*nelyU+tkyb:(uC+1)*nelyU];
    cvxb=[cvxb coordX(uC)-1:coordX(uC)+tkxb];
end

lcDof = unique(cMat(unique([reshape(elNrs(: , [1:tkzb,end-tkzb+1:end],end).',1,[]), ...
    reshape(elNrs([1:tkyb,cty,cty1,end-tkyb+1:end] , :,end).',1,[])]),[4,7,16,19]));%node numbers on loading is 2,3,6,7 x direction dofs only
lcDofix = unique(cMat(unique([reshape(elNrs(: , [1:tkzb,end-tkzb+1:end],end).',1,[]), ...
    reshape(elNrs([1:tkyb,cty,cty1,end-tkyb+1:end] , :,end).',1,[])]),[5,6,8,9,17,18,20,21]));%node numbers on loading is 2,3,6,7 y,z direction dofs only
fixed = unique(cMat(unique([reshape(elNrs(: ,[1:tkzb,end-tkzb+1:end],1).',1,[]),...     %node numbers on BC is 1,4,5,8 all dofs
     reshape(elNrs([1:tkyb,cty,cty1,end-tkyb+1:end] , :,1).',1,[])]),[1:3,10:15,22:24]));
edgedof=unique(cMat(unique([elNrs(: , [1,end],end); (elNrs([1,end], : ,end))']),[4,7,16,19]));
modF = 1/(length(lcDof)-1);                                   % modulus of the force density

%Passive regions
%columns 
pa1 = elNrs([1,[coordY(:)]',[coordY(:)]'+1,nely], [1,end,tkz,end-tkz+1] ,:);% outer edge
pa2 = elNrs([tky,[coordY(:)]'-tky+1,[coordY(:)]'+tky,nely-tky+1], [1,end] ,:);% outer along y
% beams
pa3 = elNrs(:, [1:tkzb,end-tkzb+1:end], [1:tkxb,ctx,ctx1,nelx-tkxb+1:nelx]);%top bottom beams 4mm thickness
pa4 = elNrs([1:tkyb,cty,cty1,nely-tkyb+1:nely],:,[1:tkxb,ctx,ctx1,nelx-tkxb+1:nelx]);%front back beams 4mm thickness
%Void regions
%center and faces 
pa5 = elNrs(cvy, tkz+1:end-tkz, :);% interior void
pa6 = elNrs(cvyf, [1:tkz end-tkz:end], cvxf);%face voids
pa7 = elNrs([1:tky cvy1 end-tky:end], tkz+1:end-tkz, cvxf);%face voids
% if beams to be square
pa8 = elNrs(cvyb, tkzb+1:end, [1:tkxb cvxb nelx-tkxb+1:nelx]);% face void
[pasS,pasV] = deal(unique([pa1(:);pa2(:);pa3(:);pa4(:)]),unique([... 
    pa5(:);pa6(:);pa7(:);pa8(:)]));

