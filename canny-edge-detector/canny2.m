function [ outImage ] = canny( inImage )
%canny edge detector
    PR = 2; %plot rows
    PC = 4; %plot cols
    stevec = 1;

    im=imread(inImage);
    %gaus=CalcGauss(1.4);
    gaus = fspecial('gaussian', 5, 1.4);
    G=imfilter(im,gaus);        
    
    subtightplot(PR,PC,stevec);
    imshow(im);
    stevec = stevec + 1;
    subtightplot(PR,PC,stevec);
    imshow(G);
    stevec = stevec + 1;
    
    Gy=fspecial('sobel');
    Gx=-Gy';
    
    %Gx=[-1 0 1;-1 0 1;-1 0 1];
    %Gy=[1 1 1;0 0 0;-1 -1 -1];
    Gx = [1 2 1; 0 0 0; -1 -2 -1];
    GY = [1 0 -1; 2 0 -2; 1 0 -1];
    
    GimX=imfilter(G,Gx);
    GimY=imfilter(G,Gy);
    
    mag = sqrt(double(GimX.^2+GimY.^2));
    [mag, angles] = imgradient(G,'prewitt');
    %mag = GimX+GimY;
    
    subtightplot(PR,PC,stevec);
    imshow(uint8(mag));
    stevec = stevec + 1;
    
    tl=uint8(10);
    th=uint8(20);
    
    %nl=zeros(size(G,1),size(G,2));
    %nh=zeros(size(G,1),size(G,2));
    %visited=nh;
    
    padX = 10;
    padY = 10;
    GimX = padarray(GimX, [padX padY], 0);
    GimY = padarray(GimY, [padX padY], 0);
    mag = padarray(mag, [padX padY], 0);
    angles = padarray(angles, [padX padY], 0);
    
    %angles = atan(double(abs(GimY))./double(abs(GimX))) .* 180 ./ pi;
    angles = atan2(double(GimY), double(GimX)) .*180 ./ pi;
    angles = round(angles);
    newAngles = angles;
    
    mag1 = zeros(size(mag,1));
    %mag1=mag;
    magOld = mag;
    
    radius = 1;
    
    %non-maximum suppresion
    for i=1+padX:size(mag,1)-padX
        for j=1+padY:size(mag,2)-padY
            quo = angles(i,j) / 45;
            rem = mod(angles(i,j), 45);
            if rem >= 23
                quo = ceil(quo);
            else
                quo = floor(quo);
            end
            ang = 45 * quo;
            
            ang = mod(ang,180);
            newAngles(i,j)=ang;
            
            tmp = submatrix(mag,i,j,radius);    %dobimo podmatriko; okrog piksla je radius pikslov
            
            repressed = uint8(0);
            
            if ang == 0
                % vodoravno
                tmp = tmp(radius+1,:);  %sredinska vrstica
                [maxEl, maxIdx] = max(tmp);
                
                if maxIdx == radius+1   %trenutni piksel je najvecji                    
                    mag1(i,j) = mag(i,j);
                else                    
                    mag1(i,j) = repressed;
                end
            elseif ang == 45
                %prva diagonala                
                tmp = fliplr(tmp);  %zrcalimo matriko da dobimo pravilno diagonalo
                tmp = diag(tmp);                
                [maxEl, maxIdx] = max(tmp);
                
                if maxIdx == radius+1   %trenutni piksel je najvecji                    
                    mag1(i,j) = mag(i,j);
                else                    
                    mag1(i,j) = repressed;
                end
            elseif ang == 90
                %navpicno                
                tmp = tmp(:,radius+1);  %sredinski navpicni stolpec
                [maxEl, maxIdx] = max(tmp);
                
                if maxIdx == radius+1   %trenutni piksel je najvecji                    
                    mag1(i,j) = mag(i,j);
                else                    
                    mag1(i,j) = repressed;
                end
            elseif ang == 135
                %druga diagonala                
                tmp = diag(tmp);    % tu ni potrebno zrcaliti
                [maxEl, maxIdx] = max(tmp);
                
                if maxIdx == radius+1   %trenutni piksel je najvecji                    
                    mag1(i,j) = mag(i,j);
                else                    
                    mag1(i,j) = repressed;
                end
            end
        end
    end
    
    [mag1, krneki] = nonmaxsup(mag, newAngles, 1.2);
    
    subtightplot(PR,PC,4);
    imshow(uint8(mag1));
    
    strong = uint8(255);
    weak = uint8(100);
    
    edges = zeros(size(mag,1),size(mag,2));    
    
    for i=1+padX:size(mag1,1)-padX
        for j=1+padY:size(mag1,2)-padY
            cur = mag1(i,j);
            if cur < tl
                mag1(i,j)=0;
            elseif cur < th
                %weak
                mag1(i,j)=weak;
            else
                %strong
                mag1(i,j)=strong;
            end
        end
    end    
    
    subtightplot(PR,PC,5);
    imshow(uint8(mag1));
    
    
%     for i=1+padX:size(mag1,1)-padX
%         for j=1+padY:size(mag1,2)-padY
%             cur = mag1(i,j);
%             if cur == strong
%                 sub = submatrix(mag1,i,j,1);    %ce je v okolici weak roba, strong
%                 if ismember(weak,sub)
%                     %mag1(i,j)=strong;   %weak postane strong
%                     for x=i-1:i+1
%                         for y=j-1:j+1
%                             if mag1(x,y)==weak
%                                 mag1(x,y)=strong;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end    

%     for i=1+padX:size(mag1,1)-padX
%         for j=1+padY:size(mag1,2)-padY
%             cur = mag1(i,j);
%             if cur == weak
%                 sub = submatrix(mag1,i,j,1);    %ce je v okolici weak roba, strong                               
%                 if ismember(strong,sub)
%                     mag1(i,j)=strong;
%                 else                    
%                     mag1(i,j)=0;
%                 end               
%             end
%         end
%     end 

    for i=1+padX:size(mag1,1)-padX
        for j=1+padY:size(mag1,2)-padY
            cur = mag1(i,j);
            if cur == strong
                sub = submatrix(mag1,i,j,1);    %ce je v okolici weak roba, strong                               
                if ismember(weak,sub)
                    mag1=replace(mag1,i,j,weak,strong);
                end                    
            end
        end
    end 
    
    %mag1 = hysthresh(mag,strong,weak);
    %subtightplot(PR,PC,6);
    %imshow(mag1);
    
    subtightplot(PR,PC,6);
    imshow(uint8(mag1));
    
    
    %the original image, 
    %the image after detecting edges and
    %the final image after edge linking
%     subplot(1,4,1);
%     imshow(magOld);
%     subplot(1,4,2);
%     imshow(mag);
%     subplot(1,4,3)
%     imshow(uint8(mag1));
%     subplot(1,4,4)
%     imshow(uint8(magOld)-uint8(mag1));
    bw=edge(im,'canny');
    
    %subtightplot(PR,PC,7);
    %imshow(bw);
    
    %M=sqrt(power(gx,2)+power(gy,2))
    %d(x,y) = arctan(gy/gx)
    
    %gx = -1, -1,-1;0,0,0;1,1,1
    %gy = -1,0,1;-1,0,1;-1,0,1
    outImage=mag1;
    print(mag1,'-dpng',strcat(inImage,'_edge','.png')) % print figure to a file

end

function [mat] = replace(mag,i,j,w,s)    
    for x=i-1:i+1
        for y=j-1:j+1
            if mag(x,y) == w
                replace(mag,x,y,w,s);
                mag(x,y)=s;                
            end
        end
    end
    mat=mag;                
end

