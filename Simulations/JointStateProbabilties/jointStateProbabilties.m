%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script launches simulations of joint state probabilities for a two
% leaf "tree". The two leafs are each in different location. Coalescence can
% only happen once the two leafs are in the same location. The joint 
% probabilities used are then conditioned on the coalescence history, i.e.
% the probability of the the two leafs being jointly in certain locations
% while they did not coalesce: P(l1=a,l2=b,not coalesced)
% These estimates are then compared to estimates where the joint states are
% simulated using the approach described in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% defines the per leaf probability of changing location
migration_rate = 0.1;

% defines the pairwise coalescence rate of two leafs in the same
% location
coalesence_rate = 1;

% The initial probabilities of the leafs being in a location
% to read as [l1=1,l1=2,l2=1,l2=2]
initial_states = [1 0 0 1];

% Get all possible combinations of two lineages with two different 
% locations
nrlins_tmp = (dec2bin(0:2^2-1));
nrlins = zeros(size(nrlins_tmp));
for a = 1 : size(nrlins,1)
    for b = 1 : size(nrlins,2)
        nrlins(a,b) = str2double(nrlins_tmp(a,b));
    end
end
clear a b

% find the initial configuration of locations of the two leafs
initial_configuration = zeros(1,2);
for j = 1 : 2
    tmp = initial_states(j:2:end);
    initial_configuration(tmp==1) = j-1;
end
clear j tmp

% find which of the possible configuration is the initial one, i.e. find 
% the index
start = zeros(size(nrlins,1),1);
for j = 1 : size(nrlins,1)
    if isequal(nrlins(j,:),initial_configuration)
        start(j) = 1;
        break;
    end
end   
clear j

% find for which configurations coalescence is possible
sums = zeros(size(nrlins,1),1);
for i = 1 : length(sums)
     [a,~]=hist(nrlins(i,:),2);
     a = a-1;
     a(a<0)=[];
     sums(i) = sum(a);
end
clear i a

% find out how the different joint state probabilities are connected by
% only one migration event. The probability of more than one migration
% event at a time is assumed to be sufficiently small the neglect
connectivity = zeros(length(sums),length(sums));
for a = 1 : size(connectivity,1)
    for b = 1 : size(connectivity,2)
        tmp = sum(nrlins_tmp(a,:) ~= nrlins_tmp(b,:) );
        if tmp == 1
            connectivity(a,b)=1;
        end
    end
end
clear a b tmp nrlins_tmp

% fill the diagonal elements with the column sums (matrix is symmetric
% since we assume symmetric migration)
for a = 1 : size(connectivity,1)
    connectivity(a,a)=-sum(connectivity(a,:));
end
clear a

% multiply the connectivity matrix with the migration rate already here
% such that it is not necessary late
connectivity = connectivity.*migration_rate;

% set the time span for which simulations are performed. The odes are then
% solved until the last entry of t
t=logspace(-1,1.6,50);

% simulate the joint state probabilities conditioned on survival for all
% the different time points in t 
jointStateProbs = StochSim(migration_rate,coalesence_rate,t);

% do the simulation of the joint state proabilities via the odes
[time,y] = ode45(@(time,y) conditionalODE(time,y,coalesence_rate,sums,connectivity),...
    [0 max(t)], start);

% to avoid super small values for plotting just throw the values before
% t(1) out
y(time<min(t),:)=[];
time(time<min(t))=[];

% define the colors used (standard colors)
colors = [0    0.4470    0.7410;
          0.8500    0.3250    0.0980;
          0.9290    0.6940    0.1250;
          0.4940    0.1840    0.5560];

plot(t,squeeze(jointStateProbs(1,1,:)),'o','MarkerSize',6,...
    'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:)); hold on
plot(t,squeeze(jointStateProbs(2,1,:)),'o','MarkerSize',6,...
    'MarkerEdgeColor',colors(4,:),'MarkerFaceColor',colors(2,:)); hold on
plot(t,squeeze(jointStateProbs(1,2,:)),'o','MarkerSize',6,...
    'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:)); hold on
plot(t,squeeze(jointStateProbs(2,2,:)),'o','MarkerSize',6,...
    'MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(4,:)); hold on

plot(time,y(:,1),'-','LineWidth',2,'Color',colors(1,:));hold on
plot(time,y(:,2),'--','LineWidth',2,'Color',colors(2,:));hold on
plot(time,y(:,3),'--','LineWidth',2,'Color',colors(3,:));hold on
plot(time,y(:,4),'--','LineWidth',2,'Color',colors(4,:));hold on

legend('P(l1=1,l2=1)(t) simulated','P(l1=2,l2=1)(t) simulated',...
        'P(l1=1,l2=2)(t) simulated','P(l1=2,l2=2)(t) simulated',...
        'P(l1=1,l2=1)(t) analytical','P(l1=2,l2=1)(t) analytical',...
        'P(l1=1,l2=2)(t) analytical','P(l1=2,l2=2)(t) analytical');


ylabel('joint probability of beeing in a configuration and not having coalesced')
xlabel('time')
set(gca,'YScale','log')
set(gca,'XScale','log')

