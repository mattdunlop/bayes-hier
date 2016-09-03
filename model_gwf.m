function data = model_gwf(U)

    % Model parameters
    Dx=6;               % Domain = [0,Dx]*[0,Dy]
    Dy=6;
    J=64;               % Number of measurements
    Nx=size(U,1);       % Number of x-grid points (read from coefficient)
    Ny=Nx;              % Number of y-grid points
    
    % Construct the recharge term
    Recharge = zeros(Nx,Ny);
    Recharge(ceil(4*Nx/6):Nx,:) = 137;
    Recharge(floor(5*Nx/6):Nx,:) = 274;
    Recharge = reshape(Recharge,1,Nx*Ny)/(Nx*Ny);
    
    
    [measure_x,measure_y]=meshgrid(Nx/(sqrt(J)+1):Nx/(sqrt(J)+1):Nx*(1-1/(sqrt(J)+1)),Ny/(sqrt(J)+1):Ny/(sqrt(J)+1):Ny*(1-1/(sqrt(J)+1)));
    measure_x=round(reshape(measure_x,1,J));
    measure_y=round(reshape(measure_y,1,J));
    Well_Location.pressure.x=measure_x;
    Well_Location.pressure.y=measure_y;
    
    Boundary.BN_x=zeros(Ny,2);
    Boundary.BN_y=zeros(Nx,2);
    for i=1:Ny
        Boundary.BN_x(i,2)=-Dy/Ny*500 ; %Dy=1;
    end
    for i=1:Nx
        Boundary.PD(i)=100;
    end




    %%%%%%%%%%%%%%%%%%Measurements%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Grid.Nx=Nx; Grid.hx=Dx/Nx;
    Grid.Ny=Ny; Grid.hy=Dy/Ny;
    Grid.V=Grid.hx*Grid.hy;
    Grid.N=Grid.Nx*Grid.Ny;
    Well_Data.Well_Location.pressure=Well_Location.pressure;
    Well_Data.Num_Press=length(Well_Data.Well_Location.pressure);
    Production.Q=Recharge;
    Production.BN_x=Boundary.BN_x;
    Production.BN_y=Boundary.BN_y;
    Production.PD=Boundary.PD;
    Model.Production=Production;
    Model.Grid=Grid;
    Model.Well_Data=Well_Data;

    
    Grid=Model.Grid;
    Nx=Grid.Nx;
    Ny=Grid.Ny;
    N=Nx*Ny;
    Model.Rock.K=reshape(U,Nx,Ny);
    Model.P=zeros(N,1);
    Model.Q=zeros(length(Model.Well_Data.Well_Location.pressure),1);


    Production=Model.Production;
    Rock=Model.Rock;

    % Compute TPFA
    Nx=Grid.Nx; Ny=Grid.Ny;  N=Nx*Ny;
    hx=Grid.hx; hy=Grid.hy; 
    L=(Rock.K);

    tx=hy/hx; TX=zeros(Nx+1,Ny);
    ty=hx/hy; TY=zeros(Nx,Ny+1);
    Average_x=0.5*(L(1:Nx-1,:)+L(2:Nx,:));
    TX(2:Nx,:)=Average_x.*tx;
    Average_y=0.5*(L(:,1:Ny-1)+L(:,2:Ny));
    TY(:,2:Ny)=Average_y.*ty;
    TY(1:Nx,1)=ty.*L(:,1)*2;

    x1=reshape(TX(1:Nx,:),N,1); x2=reshape(TX(2:Nx+1,:),N,1);
    y1=reshape(TY(:,1:Ny),N,1); y2=reshape(TY(:,2:Ny+1),N,1);

    for i=1:Nx
        f(i)=-y1(i)*Production.PD(i);
    end

    DiagVecs=[-y2,-x2,x1+x2+y1+y2,-x1,-y1];
    DiagIndx=[-Nx,-1,0,1,Nx];
    A=spdiags(DiagVecs,DiagIndx,N,N);


    for i=1:Grid.Nx
        Production.Q(i)=Production.Q(i)-f(i);
    end
    for j=1:Grid.Ny
        Production.Q(1+(j-1)*Grid.Nx)=Production.Q(1+(j-1)*Grid.Nx)-Production.BN_x(j,2);
        Production.Q(Grid.Nx+(j-1)*Grid.Nx)=Production.Q(Grid.Nx+(j-1)*Grid.Nx)-Production.BN_x(j,1);
    end
    for j=1:Grid.Nx
        Production.Q(j+(Grid.Ny-1)*Grid.Nx)=Production.Q(j+(Grid.Ny-1)*Grid.Nx)-Production.BN_y(j,2);
    end

    Model.P=A\Production.Q';

    u=reshape(Model.P,Nx,Ny);
    data=[];
    for i=1:1:length(Model.Well_Data.Well_Location.pressure.x)
        data=[data;u(Model.Well_Data.Well_Location.pressure.x(i),Model.Well_Data.Well_Location.pressure.y(i))];
    end

end
