function gril=grille(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz)

%--------------------------------------------------------------------------
% Function to create a regular 3D grid.
%--------------------------------------------------------------------------
%
% INPUT
%      xmin         : coordinate of first x grid node.
%      xmax         : coordinate of last x grid node.
%      dx           : grid mesh along x coordinates.
%      ymin         : coordinate of first y grid node.
%      ymax         : coordinate of last y grid node.
%      dy           : grid mesh along y coordinates.
%      zmin         : coordinate of first z grid node.
%      zmax         : coordinate of last z grid node.
%      dz           : grid mesh along z coordinates.
%--------------------------------------------------------------------------

x=[xmin:dx:xmax]';
y=[ymin:dy:ymax]';
z=[zmin:dz:zmax]';
nx=length(x);
ny=length(y);
nz=length(z);
gril=[kron(ones(ny*nz,1),x), kron(kron(ones(nz,1),y),ones(nx,1)), kron(z,ones(nx*ny,1))];
end