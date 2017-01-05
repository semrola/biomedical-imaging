function [ output ] = beat( record )
    try
        fpsN=rdann(record, 'atr',[],[],[],'N');   % manual R points - normal            
    catch exception
        fprintf('Ni N-jev\n');
        fpsN=[];        
    end
    
    try
        fpsV=rdann(record, 'atr',[],[],[],'V');   % manual R points - ventricular
    catch exception
        fprintf('Ni V-jev\n');
        fpsV=[];
    end
    
    fps=sort([fpsN' fpsV']);
    
    all=load(strcat(record,'m.mat'));
    x1 = all.val(1,:);
    x2 = all.val(2,:);
    
    first=x1;
    x1=x1+x2;
    
x1=HPFilter(x1,2,1/200);
%x1=t3+s3;
    
    %konstante
    M=9; %averaging window - probaj 12 oz. 13
    N=3;    % primeri ki jih ne upostevamo (min in max), probaj 3 oz. 4
    K=10; % 30 je za 200Hz sampling rate, 54 recimo za 360
    
    f=zeros(M,1);  % vektor utezi oz. koeficienti
    
    % trimmed moving average process - TMAF
    for k=N+1:M-N
        f(k,1)=1/(M-2*N);
    end
        
    stZaQrs=1;
    matchingNs=0;
    matchingVs=0;    
    
    Fs=200;
    T=1/Fs;
    %FP - 70 ms to FP + 140 ms
    ms140=round(0.14/T);
    ms70=round(0.07/T);
    
    output=[];        
    z=zeros(1,length(x1));
    
    
    %damo cez filtre
    for i=K:length(x1)
    %for i=K:pts        
        
        % NHPF ----> nonlinear high-pass filter
        vrednosti=x1(i-M+1:i)';  % vrednosti v oknu
        tmaf=f'*sort(vrednosti);    % koef. mnozimo s sortiranimi vrednostmi
        y(i)=x1(i-((M+1)/2))-tmaf;
        % ------>
        
        % NLPF ----> nonlinear low-pass filter
        z(i)=0;
        for k=0:K-1           
            z(i)=z(i)+power(y(i-k),2);            
        end
        % ------>        
                                       
    end
    
    thbs=[];
    thb=mean(z(50:1000));
    thb1=thb;
    beta=0.05;   %forgetting factor - med 0 in 1
    
    %dolocimo treshold
    for j=1:length(z)
        thbs=[thbs thb];    %belezimo tresholde
        lw=j-ms70;
        hi=j+ms140;
        if(lw <= 0); lw=1; end
        if(hi > length(z)); hi=length(z); end        
    
        peak=max(z(lw:hi));
        
        if(z(j) >= thb)
            thb=beta*(0.4*peak)+(1-beta)*thb;
        end
    end
    
    %gremo cez qrs-je    
    stN=[];
    stV=[];
    thbs1=[];   %za izracun thbjev v 2. zanki
    thbs2=zeros(1,length(z));   %razsirjen vektor cez celoten posnetek
    
    for i=1:length(fps)
        thbs1=[thbs1 thb1];
        idx=fps(i); %dobimo index QRSja (v x1 funkciji)
        
        thbs2(idx)=thb1;
        
        lw=idx-ms70;
        hi=idx+ms140;
        
        if(lw <= 0) lw=1; end
        if(hi > length(z)) hi=length(z); end
        
        [qrsZ,qrsZIdx]=max(z(lw:hi));   %dobimo QRS na z funkciji, gledamo okrog zaradi delaya
        qrsIdx=qrsZIdx+idx-1;
        
        if(qrsIdx > length(thbs)) qrsIdx=length(thbs); end
        thb=thbs(qrsIdx);                       
        
        round([qrsZ thb qrsIdx qrsZIdx lw idx hi]);
        
        if(qrsZ >= thb) %normal beat
            output(i)='N';
            stN=[stN idx];
            if(ismember(idx,fpsN))
                matchingNs=matchingNs+1;
            end
        else  %abnormal beat
            output(i)='V';
            stV=[stV idx];
            if(ismember(idx,fpsV))   %samo za matching
                matchingVs=matchingVs+1;
            end
        end
    end
    
    pts=304000;
    %stV(1:20);
    setdiff(fpsV,stV)
    tpts=pts;
    K=300000;
    if pts > length(thbs) tpts=length(thbs); end
    if pts > length(thbs1) tpts=length(thbs1); end
    if pts > length(thbs2) tpts=length(thbs2); end
    subplot(2,1,1);
    plot(x1(K:pts),'r');
    %subplot(4,1,2);
    %plot(first(1:pts));
    %hold on;
    subplot(2,1,2);
    plot(thbs(K:pts),'b');
    hold on;
    plot(z(K:pts),'r');
    %subplot(4,1,4);
    %plot(x2(1:pts));
    %hold on;
    %plot(x1(fps(1:1000)),'g');
    
    
    %-----PISANJE CLS---------
    
    fprintf('%3s %6s | %6s| %5s | %4s\n','','My','Correct','Match','Diff');
    %fprintf('‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n');
    fprintf('Ns: %6d | %6d | %5d | %4d\n',sum(output=='N'),size(fpsN,1),matchingNs,size(fpsN,1)-matchingNs);
    fprintf('Vs: %6d | %6d | %5d | %4d\n',sum(output=='V'),size(fpsV,1),matchingVs,size(fpsV,1)-matchingVs);
    %fprintf('SUM: %d Correct:%d\n',numel(output),size(fpsV,1)+size(fpsN,1));
    
    asciName = sprintf('%s.cls',record);
    fid = fopen(asciName, 'wt');

    for i=1:size(output,2)
        fprintf(fid,'0:00:00.00 %d %s 0 0 0\n', fps(i), output(i) );
    end
    fclose(fid);
    
    %report=bxb(record,'atr','cls','bxbReport.txt');    
end

