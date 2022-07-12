function [threshold] = Otsu_Hist(inHistogram,lowerBound, upperBound, masked)
% This function receives histogram as input. It computes a threshold 
% considering only the values within the range defined by lowerBound and
% upperBound
histN = upperBound - lowerBound + 1;
histogram = zeros(1, histN);
for i=1:histN
    histogram(i) = inHistogram(i+lowerBound-1);
end
% for i=1:lowerBound-1
%     histogram(i) = 0;
% end
% 
% for i = upperBound+1:length(histogram)
%     histogram(i) =0;
% end
% histN = upperBound;
threshold=1;
w0 = 0.0;
w1 = 0.0;
m0 = 0.0;
m1 = 0.0;
N = 0;
p = 0;
sum = 0.0;
mean = 0.0;
var_bet_class = 0.0;
var_max = 0.0;
mu_k = 0.0;

if masked ==1
    histogram(1) = 0;
end

for i=1:histN
    N=N+histogram(i);
    sum=sum+histogram(i)*i;
end

mean=sum/N;

for i=1:histN
    p=histogram(i)/N;
    
    %cumulative for class 1 and class 2
    w0=w0+p;
    w1=1-w0;
    
    %mean for class 1 and class 2
    mu_k=mu_k+i*p;
    
    if(w0==0)
        m0=0;
    else
        m0=(mu_k/w0);
    end
    
    if(1-w0==0)
        m1=0;
    else
        m1=(mean-mu_k)/(1-w0);
    end
    
    var_bet_class = w0*(m0-mean)*(m0-mean)+w1*(m1-mean)*(m1-mean);
    
    if var_bet_class>var_max
        var_max=var_bet_class;
        threshold=i;
    end
end
threshold = threshold + lowerBound - 1;
end

