function [action,E_utility_max] = best_respond(...
    nr_agents,agent,action_space,action_dimension,eed,distance)
% 'starting'
E_utility_max = 0;
%% Target assignment
% NOTE: Assumes all agents start from origin
% Compute expected distance from origin for each target
epsilon = 0; % importance of accelaration
action = [];
E_utility_agent = 0;

eeds = zeros(action_dimension,nr_agents-1);
for aa = 1:nr_agents
    if aa ~=agent
        eeds(:,aa) = eed(1,:,aa);
    end
end

for ii = 1:action_dimension
    act = action_space(ii);
    sums = 0;
    others_action = action_space;        
    E_utility_agent = prod(1-eeds(ii,:))/((distance(ii)^2));   
    if E_utility_max < E_utility_agent
        E_utility_max = E_utility_agent;
        action = act;
    elseif E_utility_max == E_utility_agent
        action = [action act];
    end
end
E_utility_max;
if E_utility_max == 0
    action = action_space(randi(action_dimension)); % select randomly if no action is available
else
[dadidadi dim] = size(action);
action = action(randi(dim)); %select randomly among equally good actions
end

end
% estimates_empirical_distribution(:,:,agent),belief_theta_mean(agent)

% do asymmetric weights - still potential function randomly generate
% do example in Brian Swenson article
% add learning theta






