function [ fltd ] = HPFilter( sig, Fc, T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    c1=1/(1+tan(Fc*pi*T));
    c2=(1-tan(Fc*pi*T))/(1+tan(Fc*pi*T));
    sigLen=size(sig, 2);
    fltd=zeros(1,sigLen);
    
    for idx=2:sigLen
        fltd(idx)=c2*fltd(idx-1)+c1*(sig(idx)-sig(idx-1));
    end
end
    
%odstranjevanje drifta + popravljanje faze    
%load('s20501m.mat')
%fltd= HPFilter(val(2,:), 2.2, 1/250)
%fliped=fliplr(fltd);
%fliped=HPFilter(fliped, 2.2, 1/250)
%ffliped=fliplr(fliped);