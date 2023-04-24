function [x_next, K, P_next] = estimate(P,H,R,x,y)

    K = P*H' / (H*P*H'+R);

    x_next = x + K * (y- H*x);


    row  = size(K,1);
    col = size(H,2);

    P_next = (eye(row, col) - K*H) * P';

%     P_next = P - K*H*P - P*H'*K'+ K * (H*P*H' + R)*K';
end