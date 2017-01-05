function [ rtrn ] = kueres2(record)
%"A QRS Complex Detection Algorithm Using Electrocardiogram Leads"

s = load(strcat(record,'m.mat'));
x1 = s.val(1,:);
x2 = s.val(2,:);
x1=HPFilter(x1, 2.2, 1/360);
x2=HPFilter(x2, 2.2, 1/360);

slen=size(x1,2);
rtrn=[];

t=cputime();

s1=zeros(1,slen);
t1=zeros(1,slen);
% first filter - low pass filter
for i=3:slen
    s1(i) = (x1(i-2) + 2*x1(i-1) + x1(i))/4;
    t1(i) = (x2(i-2) + 2*x2(i-1) + x2(i))/4;
end

s2=zeros(1,slen);
t2=zeros(1,slen);
% second filter - mains filter
for i=3:slen
    s2(i) = s1(i-2) - 2*cos((60*pi)/125) * s1(i-1) + s1(i);
    t2(i) = t1(i-2) - 2*cos((60*pi)/125) * t1(i-1) + t1(i);
end

s3=zeros(1,slen);
t3=zeros(1,slen);
% third filter - derivative filter
for i=7:slen
    s3(i) = s2(i) - s2(i-6);
    t3(i) = t2(i) - t2(i-6);
end

s3=abs(s3);
t3=abs(t3);

avgS=mean(s3);
avgT=mean(t3);

z=s3+t3;

pts=1500;

% end of signal conditioner

% 650000 vzorcev, 30min 5s = 1805s => vzorec vsake 0.00278s (360,111Hz)

% main
btMin=100;
dtMin=100;
btMax=1200;
dtMax=1200;
et=100;
Fs=360;
T=1/Fs;
step=round(0.18/T);    %180ms
ms200=round(0.2/T);    %200ms
ms16=round(0.016/T);   %16ms
ms1000=360;  %1000ms = 1s
ms500=ms1000/2;  %500ms

[najdi,najdiIdx]=max(z(1:5000)); %najdi peak v prvih n tockah    
btMin=najdi*0.405;
dtMin=btMin;

btMax=najdi*0.71;
dtMax=btMax;

bt=btMin;
dt=dtMin;

lastQRSes=[najdiIdx];   %indexi zadnjih 8 qrsjev
st=1; %stevilo rezultatov
it=0; %stevilo dejanskih iteracij
i=12; %prvih 11 primerov spustimo
stevec=0;
%slen=pts;
while i<=slen
    qrs=false;
    noise=false;
    it=it+1;
    stop=i+step;
    if stop > slen  %da ne gremo cez mejo
        stop=slen;
    end
    
    % 200ms after detection of QRS
    if((lastQRSes(end) + ms200) > i)
        et=0.2*et;
    end
    
    % no detection after 1s
    if((lastQRSes(end) + ms1000) > i)
        et=0.5*et;
    end
    
    crossings=0;
    
    if (z(i) >= dt)   % crossing zaznan
        peaks=[];
        peakIdxs=[];
        lastUpCrossingIndex=i;
        crossings=crossings+1;
        iscem=0;  %iscem < dt - potrebno za crossing ki gre pod mejo
        x=i+1;
        while x <= stop  % 180ms
            trenutni=z(x);
            prejsnji=z(x-1);
            if(trenutni < prejsnji && iscem==0 && trenutni >= dt)
                peaks=[peaks prejsnji];
                peakIdxs=[peakIdxs x-1];
                iscem=1;
            end
            
            if(trenutni < dt && iscem==1)
                lastDownCrossingIndex=x;
            end
            
            if(trenutni > prejsnji && iscem==1)
                iscem=0;                
            end
            
            if(z(x) > dt)
               stop=x+step;  %ce je crossing gledamo 180ms naprej
               if(stop > slen) stop=slen; end
            end
            x=x+1;
        end
        
        %-------------------------------------------------------------
        
        crossings=size(peaks, 2);        
        [peak, maxPeakIndex]=max(peaks);
        peakIndex=peakIdxs(maxPeakIndex);
        
        round([i, bt, dt, crossings]);
        
        %adjust bt and dt
        
        %noise
        if (crossings > 4)
            bt=1.5*bt;
            dt=max(0.5*peak,bt);
            noise=true;
        end
        
        %qrs complex
        if ((crossings>=2) && (crossings<=4))
            rtrn(st)=peakIndex;%lastUpCrossingIndex;
            st=st+1;
            bt = (0.75*bt + 0.25 * peak);
            dt = max(0.5 * peak, bt);
            qrs=true;
        end
        
        %secondary detector
        if (crossings == 1)
            %rtrn(st)=peakIndex;
            %st=st+1;
            stevec=stevec+1;
            %lahko tudi s symsum
            vsota=0;
            for idx=0:19
                if (i-idx) > 0
                    vsota=vsota+power(z(i),2)*(i-idx);
                end
            end
            energy=vsota/(idx+1);
            % energy > et
            % RR > 200ms
            % width between 16ms & 500ms
            % amplitude between 16% & 600% of mean of last 8 QRSes            
            width=lastDownCrossingIndex - lastUpCrossingIndex;
            ampAvg=mean(z(lastQRSes));
            if ampAvg > peak ampAvg=ampAvg*1.25; end
            if ampAvg < peak ampAvg=ampAvg*0.75; end
            
            ampLow=0.1*ampAvg;
            ampHi=6*ampAvg;            
                            
            if(energy >= et && (not(isempty(lastQRSes)) && peakIndex-lastQRSes(end) > ms200) ...
                    && (width >= ms16 && width <= ms500) ...
                    && (peak >= ampLow && peak <= ampHi))
                qrs=true;
                rtrn(st)=peakIndex;%lastUpCrossingIndex;
                st=st+1;               
            end
        end
        
        if(qrs)
            et=4*(0.75*et + 0.5*peak);
            %belezenje
            if (size(lastQRSes,2) < 8)
                lastQRSes=[lastQRSes peakIndex];
            else
                lastQRSes=[ lastQRSes(2:end) peakIndex ]; %dodamo nov qrs na konec
            end
        end
        
        % first critical time - no QRS complex detected
        if(not(qrs))
            bt= 0.5*bt;
            dt= bt;
        end
        
        %second critical time
        if(qrs || noise)
            dt = 0.75*dt;
            dt = max(dt, 0.5*bt);
        end
        
        i=stop+1;
        
    else  % ni crossinga
        i=i+1;
    end        
    
    if(dt < dtMin)
        dt=dtMin;
    end
    if(bt < btMin)
        bt=btMin;
    end
    if(et < btMin)
        et=btMin;
    end
    if(dt > dtMax)
        dt=dtMax;
    end
    if(bt > btMax)
        bt=btMax;
    end
    if(et > btMax)
        et=btMax;
    end
    
    crs(it)=crossings;
    dts(it)=dt;
    bts(it)=bt;
    
end

fprintf('Record: %s time: %f\n', record, (cputime()-t));
asciName = sprintf('%s.det',record);
fid = fopen(asciName, 'wt');

for i=1:size(rtrn,2)
    fprintf(fid,'0:00:00.00 %d N 0 0 0\n', rtrn(1,i) );
end
fclose(fid);

end