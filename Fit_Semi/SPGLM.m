clear all
addpath 'data_ch211'
load('train_data')
load('len')
% load data
X_all=train_sort(1:end,1);
Y_all=train_sort(1:end,2);
% normalize the data set
% pick the top 3000 points
X=X_all(1:3000);
Y=Y_all(1:3000);
Y=(sqrt(Y)-sqrt(X))./sqrt(X);
X=sqrt(X);
% We only consider Y \in [-0.8,0.8]
ind_1=find(Y<0.8);
ind_2=find(Y>-0.8);
Y=Y(intersect(ind_1,ind_2));
X=X(intersect(ind_1,ind_2));
% consider the intercept as well
X=[X,ones(length(X),1)];

%% parameters that will not change during the iterations
unique_Y=unique(Y);
K=length(unique_Y);
N=length(Y);
% difference on Y's
Delta_vec=unique_Y(2:end)-unique_Y(1:end-1);
% B is the transformation matrix from \phi to \alpha
B=zeros(K,K);
for i = 1:K
    if i ==1
        B(1,1)=1;
        B(1,2)=-1/Delta_vec(1);
    elseif i==K;
        B(K,K-1)=-1/Delta_vec(end-1);
        B(K,K)=1;
    else
        B(i,i-1)=-1/Delta_vec(i-1);
        B(i,i)=1/Delta_vec(i-1)+1/Delta_vec(i);
        B(i,i+1)=-1/Delta_vec(i);
    end
end

%% Initialization
% initialization of beta
beta_1=sum((X(:,1)-mean(X(:,1))).*(Y-mean(Y)))/sum((X(:,1)-mean(X(:,1))).^2);
beta_2=mean(Y)-mean(X(:,1))*beta_1;
beta_coe=[beta_1,beta_2]';
theta_vec=X*beta_coe;

% initialization of phi (choose optimal phi)
% firstly we generate the slope vector, which contains (K-1) slopes, we
% order them in ascending order
slope_vec=sort(rand([1,K-1]),2,'descend');
phi_vec=ones(1,K);
for k = 2:K
    phi_vec(k)=Delta_vec(k-1)*slope_vec(k-1)+phi_vec(k-1);
end
phi_vec=phi_vec(:);
% difference on phi's
delta_vec=phi_vec(2:end)-phi_vec(1:end-1);
term_1=exp(theta_vec*unique_Y');
term_2=exp(repmat(phi_vec',N,1));
%% Optimization
% Gradient Descent
f_global_current=1;
f_global_old=0;
f_global_value=[log_concave_lik(phi_vec,Y,Delta_vec,theta_vec,term_1,term_2)];
global_count=1;
soln_beta=[beta_coe'];
soln_phi=[phi_vec'];
% we stop the optimization after f_global does not change anymore or 
% certain rounds of iterations are finished
while abs((f_global_current-f_global_old)./f_global_current)>0.000001 & global_count<500
    tic
    f_current=1;
    f_old=0;
    f_beta_value=[f_global_value(end)];
    % optimization on beta_coe given fixed phi_vec
    beta_done=0;
    beta_cnt=1;
    while abs((f_current-f_old)/f_current)>0.000001 && beta_cnt<100
        step_size=1/2^10;
        cnt=1;
        J_vec=[];
        for i=1:N
            temp=0;
            for k = 1:K-1
                temp=temp+Delta_vec(k)/(delta_vec(k)+theta_vec(i)*Delta_vec(k))*(exp(theta_vec(i)*unique_Y(k+1)+phi_vec(k+1))...
                    -exp(theta_vec(i)*unique_Y(k)+phi_vec(k)));
            end
            J_vec=[J_vec,temp];
        end
        
        J_pb1=[];
        for i=1:N
            temp=0;
            for k = 1:K-1
                temp=temp+Delta_vec(k)/(delta_vec(k)+theta_vec(i)*Delta_vec(k))*(X(i,1)*unique_Y(k+1)*exp(theta_vec(i)*unique_Y(k+1)+phi_vec(k+1))-...
                    X(i,1)*unique_Y(k)*exp(theta_vec(i)*unique_Y(k)+phi_vec(k)))-...
                    Delta_vec(k)^2*X(i,1)/(delta_vec(k)+theta_vec(i)*Delta_vec(k))^2*...
                    (exp(theta_vec(i)*unique_Y(k+1)+phi_vec(k+1))-exp(theta_vec(i)*unique_Y(k)+phi_vec(k)));
            end
            J_pb1=[J_pb1,temp];
        end
        
        J_pb0=[];
        for i=1:N
            temp=0;
            for k = 1:K-1
                temp=temp+Delta_vec(k)/(delta_vec(k)+theta_vec(i)*Delta_vec(k))*...
                    (unique_Y(k+1)*exp(theta_vec(i)*unique_Y(k+1)+phi_vec(k+1))-...
                    unique_Y(k)*exp(theta_vec(i)*unique_Y(k)+phi_vec(k)))-...
                    Delta_vec(k)^2/(delta_vec(k)+theta_vec(i)*Delta_vec(k))^2*...
                    (exp(theta_vec(i)*unique_Y(k+1)+phi_vec(k+1))-...
                    exp(theta_vec(i)*unique_Y(k)+phi_vec(k)));
            end
            J_pb0=[J_pb0,temp];
        end
        J_score_beta=[sum(X(:,1).*Y)-sum(J_pb1./J_vec),sum(Y)-sum(J_pb0./J_vec)]';
        f_old=f_beta_value(end);
        beta_new=beta_coe+step_size*J_score_beta;
        theta_vec_new=X*beta_new;
        term_1_new=exp(theta_vec_new*unique_Y');
        f_current=log_concave_lik(phi_vec,Y,Delta_vec,theta_vec_new,term_1_new,term_2);
        while isnan(f_current) | f_current<f_old 
            step_size=step_size/2;
            beta_new=beta_coe+step_size*J_score_beta;
            theta_vec_new=X*beta_new;
            term_1_new=exp(theta_vec_new*unique_Y');
            f_current=log_concave_lik(phi_vec,Y,Delta_vec,theta_vec_new,term_1_new,term_2);
            cnt=cnt+1;
            if cnt>100
                break
            end
        end
        f_beta_value=[f_beta_value,f_current];
        beta_coe=beta_new;
        theta_vec=X*beta_coe;
        term_1=exp(theta_vec*unique_Y');
        if abs((f_current-f_old)/f_current)<0.000001
            beta_done=1
        end
        beta_cnt=beta_cnt+1;
    end
    % optimization on phi_vec given fixed beta_vec
    f_current=1;
    f_old=0;
    f_phi_value=[f_beta_value(end)];
    phi_cnt=1;
    phi_done=0;
    while abs((f_current-f_old)./f_current)>0.000001 & phi_cnt<100
        step_size=1/2^10; 
        cnt=1;
        J_pmat=zeros(K,N);
        for i = 1:N
            for k = 1:K
                if k>1 & k < K
                    temp=Delta_vec(k)*exp(theta_vec(i)*unique_Y(k)+...
                        phi_vec(k))*(-delta_vec(k)-theta_vec(i)*Delta_vec(k)+...
                        exp(theta_vec(i)*Delta_vec(k)+delta_vec(k))-1)/...
                        (delta_vec(k)+theta_vec(i)*Delta_vec(k))^2+...
                        Delta_vec(k-1)*exp(theta_vec(i)*unique_Y(k)+phi_vec(k))...
                        *(delta_vec(k-1)+theta_vec(i)*Delta_vec(k-1)-1+...
                        exp(-theta_vec(i)*Delta_vec(k-1)-delta_vec(k-1)))/...
                        (delta_vec(k-1)+theta_vec(i)*Delta_vec(k-1))^2;
                elseif k==1
                    temp=Delta_vec(k)*exp(theta_vec(i)*unique_Y(k)+phi_vec(k))...
                        *(-delta_vec(k)-theta_vec(i)*Delta_vec(k)+...
                        exp(theta_vec(i)*Delta_vec(k)+delta_vec(k))-1)/...
                        (delta_vec(k)+theta_vec(i)*Delta_vec(k))^2;
                else
                    temp=Delta_vec(k-1)*exp(theta_vec(i)*unique_Y(k)+phi_vec(k))*...
                        (delta_vec(k-1)+theta_vec(i)*Delta_vec(k-1)-1+...
                        exp(-theta_vec(i)*Delta_vec(k-1)-delta_vec(k-1)))/...
                        (delta_vec(k-1)+theta_vec(i)*Delta_vec(k-1))^2;
                end
                J_pmat(k,i)=temp;
            end
        end
        J_score_phi=[];
        for k = 1:K
            J_score_phi=[J_score_phi,sum(Y==unique_Y(k))-sum(J_pmat(k,:)./J_vec)];
        end
        alpha_old=B*phi_vec;
        alpha_new=alpha_old+step_size*inv(B)*J_score_phi';
        phi_vec_new=inv(B)*alpha_new;
        term_2_new=exp(repmat(phi_vec_new',N,1));
        f_old=f_phi_value(end);
        f_current=log_concave_lik(phi_vec_new,Y,Delta_vec,theta_vec,term_1,term_2_new);
        while isnan(f_current) | f_old>f_current 
            step_size=step_size/2;
            alpha_new=alpha_old+step_size*inv(B)*J_score_phi';
            phi_vec_new=inv(B)*alpha_new;
            term_2_new=exp(repmat(phi_vec_new',N,1));
            f_current=log_concave_lik(phi_vec_new,Y,Delta_vec,theta_vec,term_1,term_2_new);
            cnt=cnt+1;
            if cnt >100
                break
            end
        end
        alpha_new=B*phi_vec_new;
        % project alpha_new to the feasible set 
        ind=find(alpha_new<0);
        ind=setdiff(ind,[1,K]);
        alpha_new(ind)=0;       
        phi_vec=inv(B)*alpha_new;
        term_2=exp(repmat(phi_vec',N,1));
        f_current=log_concave_lik(phi_vec,Y,Delta_vec,theta_vec,term_1,term_2);
        f_phi_value=[f_phi_value,f_current];
        % update of delta_vec and J_vec
        delta_vec=phi_vec(2:end)-phi_vec(1:end-1);
        J_vec=[];
        for i=1:N
            temp=0;
            for k = 1:K-1
                temp=temp+Delta_vec(k)/(delta_vec(k)+theta_vec(i)*Delta_vec(k))*(exp(theta_vec(i)*unique_Y(k+1)+phi_vec(k+1))...
                    -exp(theta_vec(i)*unique_Y(k)+phi_vec(k)));
            end
            J_vec=[J_vec,temp];
        end
        phi_cnt=phi_cnt+1;
        if abs((f_current-f_old)./f_current)<0.000001
            phi_done=1
        end
    end
    time=toc;
    f_global_current=log_concave_lik(phi_vec,Y,Delta_vec,theta_vec,term_1,term_2);
    f_global_value=[f_global_value,f_global_current];
    f_global_old=f_global_value(end-1);
    dataformat = '%d-th iteration: Optimal value = %f, time = %f\n';
    dataValue = [global_count, f_global_current, time];
    fprintf(dataformat, dataValue);
    global_count=global_count+1;
    soln_beta=[soln_beta;beta_coe'];
    soln_phi=[soln_phi;phi_vec'];
end
%% modification on optimal solution
% solve theta by binary search
ind=find(f_global_value==max(f_global_value));
phi_vec=soln_phi(ind,:)';
beta_coe=soln_beta(ind,:)';
theta_vec=X*beta_coe;
delta_vec=phi_vec(2:end)-phi_vec(1:end-1);

left_sum=0;
right_sum=0;
theta=-10;
c=mean(unique_Y);
while -left_sum+5*right_sum >= 0
    left_sum=0;
    right_sum=0;
    theta=theta+1;
    for k=1:K-1
        left_sum=left_sum+(exp(theta*unique_Y(k+1)+phi_vec(k+1))-exp(theta*unique_Y(k)+phi_vec(k)))/(delta_vec(k)+theta*Delta_vec(k))*Delta_vec(k);
        term_1=((delta_vec(k)*Delta_vec(k)+theta*Delta_vec(k)^2)*unique_Y(k+1)-Delta_vec(k)^2)*exp(theta*unique_Y(k+1)+phi_vec(k+1));
        term_2=((delta_vec(k)*Delta_vec(k)+theta*Delta_vec(k)^2)*unique_Y(k)-Delta_vec(k)^2)*exp(theta*unique_Y(k)+phi_vec(k));
        right_sum=right_sum+(term_1-term_2)/(delta_vec(k)+theta*Delta_vec(k))^2;
    end
end
theta_1=theta;

left_sum=0;
right_sum=0;
theta=-10;
while -left_sum+5*right_sum <= 0
    left_sum=0;
    right_sum=0;
    theta=theta+1;
    for k=1:K-1
        left_sum=left_sum+(exp(theta*unique_Y(k+1)+phi_vec(k+1))-exp(theta*unique_Y(k)+phi_vec(k)))/(delta_vec(k)+theta*Delta_vec(k))*Delta_vec(k);
        term_1=((delta_vec(k)*Delta_vec(k)+theta*Delta_vec(k)^2)*unique_Y(k+1)-Delta_vec(k)^2)*exp(theta*unique_Y(k+1)+phi_vec(k+1));
        term_2=((delta_vec(k)*Delta_vec(k)+theta*Delta_vec(k)^2)*unique_Y(k)-Delta_vec(k)^2)*exp(theta*unique_Y(k)+phi_vec(k));
        right_sum=right_sum+(term_1-term_2)/(delta_vec(k)+theta*Delta_vec(k))^2;
    end
end
theta_2=theta;

dif=inf;
theta=(theta_1+theta_2)/2;
left_sum=1;
sol_vec=[theta];
while abs(dif/left_sum)>0.0000001
    left_sum=0;
    right_sum=0;
    for k=1:K-1
        left_sum=left_sum+(exp(theta*unique_Y(k+1)+phi_vec(k+1))-exp(theta*unique_Y(k)+phi_vec(k)))/(delta_vec(k)+theta*Delta_vec(k))*Delta_vec(k);
        term_1=((delta_vec(k)*Delta_vec(k)+theta*Delta_vec(k)^2)*unique_Y(k+1)-Delta_vec(k)^2)*exp(theta*unique_Y(k+1)+phi_vec(k+1));
        term_2=((delta_vec(k)*Delta_vec(k)+theta*Delta_vec(k)^2)*unique_Y(k)-Delta_vec(k)^2)*exp(theta*unique_Y(k)+phi_vec(k));
        right_sum=right_sum+(term_1-term_2)/(delta_vec(k)+theta*Delta_vec(k))^2;
    end
    dif=5*right_sum-left_sum;
    if dif < 0
        theta_1=theta;
        theta=(theta_1+theta_2)/2;
    else
        theta_2=theta;
        theta=(theta_1+theta_2)/2;
    end
    sol_vec=[sol_vec,theta];
end


theta=sol_vec(end-1);
phi_vec=phi_vec+theta*unique_Y-log(left_sum);
delta_vec=phi_vec(2:end)-phi_vec(1:end-1);
% check the constraints
temp_sum=0;
for k=1:K-1
    temp_sum=temp_sum+Delta_vec(k)*(exp(phi_vec(k+1))-exp(phi_vec(k)))/delta_vec(k);
end
temp_sum

temp_sum_2=0;
for k=1:K-1
    term_1=(unique_Y(k+1)*delta_vec(k)/Delta_vec(k)-1)*exp(phi_vec(k+1));
    term_2=(unique_Y(k)*delta_vec(k)/Delta_vec(k)-1)*exp(phi_vec(k));
    temp_sum_2=temp_sum_2+(Delta_vec(k)/delta_vec(k))^2*(term_1-term_2);
end
temp_sum_2

beta_coe=beta_coe+[0,-theta]';
theta_vec=X*beta_coe;
term_1=exp(theta_vec*unique_Y');
term_2=exp(repmat(phi_vec',N,1));
opt_val=log_concave_lik(phi_vec,Y,Delta_vec,theta_vec,term_1,term_2);
beta_opt=beta_coe;
phi_opt=phi_vec;

cd ../Semi_Parameters
save('blocks','X','Y')
save('phi','soln_phi')
save('beta','soln_beta')
save('opt_val','f_global_value')
save('opt_beta','beta_opt')
save('opt_phi','phi_opt')

