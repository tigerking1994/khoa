function f = PSO(G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, G11, G12, G13, G14, KTQT, GHSVL, DCX, HSc1, HSc2, w1, w2)

LB = [G4 G10 G6 G8];         %lower bounds of variables 
UB = [G5 G11 G7 G9];      %upper bounds of variables 

% pso parameters values 
m = 4;            % number of variables 
n = KTQT;          % population size 
wmax = 0.9;       % inertia weight 
wmin = 0.4;       % inertia weight 
c1 = HSc1;           % acceleration factor 
c2 = HSc2;           % acceleration factor   
maxite = GHSVL;    % set maximum number of iteration 
    
for i=1:n         
    for j=1:m             
        x0(i,j)=round(LB(j)+rand()*(UB(j)-LB(j)));         
    end
end
x = x0;            
v = 0.1*x0;        
for i=1:n         
    f0(i,1)=ofun(x0(i,:), G1, G2, G3, G12, G13, G14, w1, w2);     
end
[fmin0,index0]=min(f0);         
pbest=x0;               % initial pbest     
gbest=x0(index0,:);     % initial gbest     

ite=1;         
tolerance=1;   
while ite<=maxite && tolerance>DCX                  
    w=wmax-(wmax-wmin)*ite/maxite; % update inertial weight           
    % pso velocity updates         
    for i=1:n             
        for j=1:m                 
            v(i,j)=w*v(i,j)+c1*rand()*(pbest(i,j)-x(i,j))...                         
                +c2*rand()*(gbest(1,j)-x(i,j));             
        end
    end

    % pso position update         
    for i=1:n             
        for j=1:m                 
            x(i,j)=x(i,j)+v(i,j);             
        end
    end

    % handling boundary violations         
    for i=1:n             
        for j=1:m                
            if x(i,j)<LB(j)                     
                x(i,j)=LB(j);                 
            elseif x(i,j)>UB(j)                     
                x(i,j)=UB(j);                 
            end
        end
    end

%     x_plt(:,:,ite) = x;
%     sz = size(x)
    plt(:,:,ite) = x;
    
    % evaluating fitness         
    for i=1:n             
        f(i,1)=ofun(x(i,:), G1, G2, G3, G12, G13, G14, w1, w2);         
    end

    % updating pbest and fitness         
    for i=1:n             
        if f(i,1)<f0(i,1)
            pbest(i,:)=x(i,:);                 
            f0(i,1)=f(i,1);             
        end
    end

    [fmin,index]=min(f0);   % finding out the best particle         
    ffmin(ite,1)=fmin;    % storing best fitness         
    ffite=ite;         % storing iteration count           

    % updating gbest and best fitness         
    if fmin<fmin0             
        gbest=pbest(index,:);             
        fmin0=fmin;         
    end

    % calculating tolerance         
    if ite>100             
        tolerance=abs(ffmin(ite-100,1)-fmin0);         
    end

%     % displaying iterative results         
%     if ite==1             
%         disp(sprintf('Iteration    Best particle    Objective fun'));         
%     end
%     disp(sprintf('%8g  %8g          %8.4f',ite,index,fmin0));             
    ite=ite+1;     
end

% pso algorithm-----------------------------------------------------end     
% gbest;          
% fvalue = 0.1441*((gbest(1))^(-0.3023))*((gbest(2))^(0.2608))*((gbest(3))^(0.5277));
% fff = fvalue;     
rgbest(1,:) = gbest;     
% disp(sprintf('--------------------------------------')); 

% pso main program------------------------------------------------------end 
% disp(sprintf('\n')); 
% disp(sprintf('*********************************************************')); 
% disp(sprintf('Final Results-----------------------------')); 
% bestfun = fff;
best_variables = rgbest(1,:);
% disp(sprintf('*********************************************************')); 
% toc   

% % PSO convergence characteristic 
% plot(ffmin(1:ffite,1),'-k'); 
% xlabel('Iteration'); 
% ylabel('Function value'); 
% title('PSO convergence characteristic') 
% %##########################################################################

% PSO plot 3D
sz = size(plt);
fig = figure('visible', 'off');
for i=1:sz(3)
    mat = plt(:,:,i);
    x = mat(:,1);
    y = mat(:,2);
    z = mat(:,3);
    plot3(x,y,z,'g.','MarkerSize',20)
    title('Patides travel path in the process')
    xlabel('Speed (x)') 
    ylabel('Feed (y)') 
    zlabel('Depth of cut (z)') 
%     axis([0 5 0 10 -10 15])
    grid on
    hold on
end
Frame = frame2im(getframe(fig));
fwidth = size(Frame, 1);
fheight = size(Frame, 2);

f = {best_variables, Frame, fwidth, fheight};
end

