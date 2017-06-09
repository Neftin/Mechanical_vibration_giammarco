
%try algorithm

clear all

vect_t = 0;

x = meshgrid(-10:10,-10:10);

y = transpose(x);
    

z = x.^2 + y.^2;

surf(x,y,z);

t = [-1,100
    ] %first guess

z_old = t(1).^2 + t(2).^2 % first eval

t_old = t;

t = [2,1] %second guess

z_new = t(1).^2 + t(2).^2% second eval

t_new = t;

i = 0

impr = 10000;

alpha = 0.03;

while (impr >= 10^(-5) )
    
    
    if (i == 0)
    
    [t_new , impr] = update_discrete_gradient(z_new, z_old, t_new, z_old, alpha);

    else  
    

    [t_new , impr] = update_discrete_gradient(z_new, z_old, t_new, z_old, alpha);
   
    end
    
    z_old = z_new
    
    z_new = t_new(1)^2  + t_new(2)^2;
    
    t_old = t_new;
    
    
    i = i+1
    
    pause;
    
    vect_t = [vect_t,t_old];
    

end




