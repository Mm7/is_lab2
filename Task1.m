%% TASK 1


%Ty|x (8  possibilities)
y_pos=[0,0,0,0,0,0,0;0,0,0,0,0,0,1;0,0,0,0,0,1,0;0,0,0,0,1,0,0;0,0,0,1,0,0,0;0,0,1,0,0,0,0;0,1,0,0,0,0,0;1,0,0,0,0,0,0];
%Tz|x (64 possibilities)
z_pos=[0,0,0,0,0,0,0;0,0,0,0,0,0,1;0,0,0,0,0,1,0;0,0,0,0,1,0,0;0,0,0,1,0,0,0;0,0,1,0,0,0,0;0,1,0,0,0,0,0;1,0,0,0,0,0,0;
      0,0,0,0,0,1,1;0,0,0,0,1,0,1;0,0,0,1,0,0,1;0,0,1,0,0,0,1;0,1,0,0,0,0,1;1,0,0,0,0,0,1;0,0,0,0,1,1,0;0,0,0,1,0,1,0;
      0,0,1,0,0,1,0;0,1,0,0,0,1,0;1,0,0,0,0,1,0;0,0,0,1,1,0,0;0,0,1,0,1,0,0;0,1,0,0,1,0,0;1,0,0,0,1,0,0;0,0,1,1,0,0,0;
      0,1,0,1,0,0,0;1,0,0,1,0,0,0;0,1,1,0,0,0,0;1,0,1,0,0,0,0;1,1,0,0,0,0,0;0,0,0,0,1,1,1;0,0,0,1,0,1,1;0,0,1,0,0,1,1;
      0,1,0,0,0,1,1;1,0,0,0,0,1,1;0,0,0,1,1,0,1;0,0,1,0,1,0,1;0,1,0,0,1,0,1;1,0,0,0,1,0,1;0,0,1,1,0,0,1;0,1,0,1,0,0,1;
      1,0,0,1,0,0,1;0,1,1,0,0,0,1;1,0,1,0,0,0,1;1,1,0,0,0,0,1;0,0,0,1,1,1,0;0,0,1,0,1,1,0;0,1,0,0,1,1,0;1,0,0,0,1,1,0;
      0,0,1,1,0,1,0;0,1,0,1,0,1,0;1,0,0,1,0,1,0;0,1,1,0,0,1,0;1,0,1,0,0,1,0;1,1,0,0,0,1,0;0,0,1,1,1,0,0;0,1,0,1,1,0,0;
      1,0,0,1,1,0,0;0,1,1,0,1,0,0;1,0,1,0,1,0,0;1,1,0,0,1,0,0;0,1,1,1,0,0,0;1,0,1,1,0,0,0;1,1,0,1,0,0,0;1,1,1,0,0,0,0];

x(1:7)=uint8(input('insert the input message (in vector):'));
% Legitimate channel
t_y=randi(length(y_pos));
c_y=y_pos(t_y,:);
y=uint8(xor(x,c_y)); % y has at most 1 binary error per word

% Eavesdropper channel
t_z=randi(length(z_pos));
c_z=z_pos(t_z,:);
z=uint8(xor(x,c_z));% z has at most 3 binary error per word

disp(['x: (input message)']);
disp([x]);
disp(['y: (legitimate channel)']);
disp([y]);
disp(['z:(eavesdropper channel)']);
disp([z]);

      
%% Verification of the conditional independence and uniformity   
l_y=0;
l_z=0;
l=0;
n=10^5;
x=[1,0,0,1,0,0,0];
for i=1:n
    t_y=randi(length(y_pos));
    c_y=y_pos(t_y,:);
    y=xor(x,c_y);
    t_z=randi(length(z_pos));
    c_z=z_pos(t_z,:);
    z=xor(x,c_z);
    
    if(y==x)
        l_y=l_y+1;
    end
    if(z==x)
        l_z=l_z+1;
    end
    if(z==x & y==z)
        l=l+1;
    end
end
Py_x=l_y/n;
Pz_x=l_z/n;
Pyz_x=l/n;

%results of the verification 
txt=sprintf('Prob. of Py|x : %d,   Ideal Prob. of Py|x  (1/8)  : %d',[Py_x , 1/8]);
disp(txt); 

txt=sprintf('Prob. of Pz|x : %d,   Ideal Prob. of Pz|x  (1/64) : %d',[Pz_x , 1/64]);
disp(txt); 

txt=sprintf('Prob. of Pyz|x : %d,  Ideal Prob. of Pyz|x (1/512): %d',[Pyz_x , 1/512]);
disp(txt); 
