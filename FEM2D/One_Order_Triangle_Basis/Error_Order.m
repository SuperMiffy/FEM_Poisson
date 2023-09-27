function [orderL2,errorL2,orderH1,errorH1] = Error_Order(N1,N2)
% this is the process of L_2 order and H_1 order
% input: refer to the main function Finite_element_1d.m
% output: order_L2: the covergence order of L2 norm
%         order_H1:the convergence order of H1 norm


Number = 6; % number of refinement
errorL2 = zeros(Number,1); % L2 error
errorH1 = zeros(Number,1); % H1 error

for i=1:Number
   [errorL2(i),errorH1(i)] = run_main_triangle_2d(N1* 2^(i-1),N2* 2^(i-1));% main function
end

orderL2 = zeros(Number-1,1);%precondition for L2 convergence
orderH1 = zeros(Number-1,1);%precondition for H1 convergence

for i=1:Number-1
   orderL2(i) = log2(errorL2(i)/errorL2(i+1));
   orderH1(i) = log2(errorH1(i)/errorH1(i+1));
end
end