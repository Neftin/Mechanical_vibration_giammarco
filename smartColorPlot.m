function [ CLR ] = colorPlot(N,phase,brightness_range,mode)
% This funtion retutns out N colors based on perceiving brightness
if (N >= 9)
    display('using more than 8 colors is not recommended, colors may be similar each other. If you are an artist, GOOD CHOICE, rainbows are fucking cute')
end

%   tone

%   dividing a HUE circle in N points

step   = round(256/N);
hue   = zeros(1,N);
    for i = 0:N-1
        if   (phase + i*step >= 255) % restart from 0
        hue(i+1) = (1+phase + i*step) - 255;
        else 
        hue(i+1) = (1+phase + i*step);
        end
    end
    
%   maixmal saturation to maximize differences

sat     = ones(1,N);             %  max saturation

if mode == 'perceived'
    
    CLR     = [hue./256 ; sat ; sat]';
    
    CLR     = hsv2rgb(CLR);                  %  rgb vestors
    
    intPond = [0.3 0.59 0.11];           %  color intensities
    
    int     = CLR*intPond';               %  normalized intensities vector
    
    [~,idx] = sort(int,'ascend');        %  index for sorting
    
    hue     = hue(idx);
    
else
    if mod(N,2) ~= 0
        hue = [ hue 0 ];
    end
    
    Hmatrix = reshape(hue, [length(hue)/2 , 2] );
    HUE     = [];
    for i = 1 : length(hue)/2
        HUE = [HUE Hmatrix(i,:)];
    end
    
    hue = HUE;
    
    if mod(N,2) ~= 0
        hue = hue( 1:(length(hue)-1) );
    end
    
end

step    = brightness_range/N;

brig   = zeros(1,N);

    for i = N-1:-1:0
        brig(i+1) = (1 - (i*step));
    end

CLR = [hue./256 ; sat ; brig]';

CLR = hsv2rgb(CLR);

%   fix bright and use the sensitivity brightness 

%   best phase around 220

end

