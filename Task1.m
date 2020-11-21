%% TASK 1

x(1:7)=uint8(input('insert the input message (in vector):'));
[y,z] = wiretrap_channel(x);

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
    [y,z] = wiretrap_channel(x);
    
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
