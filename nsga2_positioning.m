function [global_solution,global_best, final_metrics, Simulation_time_NSGA2] = nsga2_positioning(CaseSelect, Pt, N_users, X, Y, h, DataRate, ...
    Users, Demand, Zcov, Zcap, Rcov, Rcap, num_generations, Tol, pop_size, alpha, beta, normal_values, plots)

    tic
    
    crossover_prob = 0.9;
    mutation_prob = 0.1;

    N_users = length(Users);

    R_limit = min(Rcov, Rcap)*10^3;
    
    N_UAV = max(Zcov, Zcap);
    
    UAV_Lb = [0 0];
    UAV_Ub = [X Y];
    
    lb = [];
    ub = [];

    % Weights for the different objectives / criteria

    switch CaseSelect
        case 1
            w = [0.5559, 0.1364, 0.0489, 0.2589];
        case 2
            w = [0.0834, 0.6259, 0.2229, 0.0678];
        case 3
            w = [0.1558, 0.0856, 0.0491, 0.7095];
    end
    
    % Simple bounds of the search domain
    
    for r = 1:1:N_UAV
    % Lower bounds
    lb = horzcat(lb,UAV_Lb);
    % Upper bounds
    ub = horzcat(ub,UAV_Ub);
    end

    % Initialize population
    population = rand(pop_size, length(lb)) .* (ub - lb) + lb;

    for fit_pop = 1:pop_size
        [fitness(fit_pop,:), Assoc(fit_pop), Mean_SINR(fit_pop), SINR_max(:,fit_pop), Mean_TP_eMBB(fit_pop)] = objective_functions(population(fit_pop,:));
    end

    global_solution = [];
    best_solutions = [];
    best_fitness_values = [];
    fmin = 0;
    global_best = [0, 0, 0, 0];
    best_Assoc = [];
    best_SINR = [];
    best_DD = [];
    final_metrics = [0, 0, 0, 0];

    gen = 1;
    
    while gen < num_generations && (- fmin < Tol || final_metrics(4) < normal_values(4))
        % Non-dominated sorting
        [fronts, ranks] = non_dominated_sorting(fitness);
        
        % Crowding distance
        crowding_distances = crowding_distance(fitness, fronts);
        
        % Selection
        mating_pool = tournament_selection(population, ranks, crowding_distances);
        
        % Crossover & Mutation
        offspring = crossover(mating_pool, crossover_prob, lb, ub);
        offspring = mutation(offspring, mutation_prob, lb, ub);
        
        % Evaluate offspring
        for fit_pop = 1:pop_size
            [offspring_fitness(fit_pop,:), Assoc_off(fit_pop), Mean_SINR_off(fit_pop), SINR_max_off(:,fit_pop), Mean_TP_off(fit_pop)] = objective_functions(offspring(fit_pop,:));
        end
        
        % Combine & Select next generation
        combined_population = [population; offspring];
        combined_fitness = [fitness; offspring_fitness];
        
        [population, fitness] = select_next_generation(combined_population, combined_fitness, pop_size);

        [best_fronts, ~] = non_dominated_sorting(fitness);

        if ~isempty(best_fronts) && ~isempty(best_fronts{1})
            valid_indices = best_fronts{1}(best_fronts{1} <= size(population, 1));
            if ~isempty(valid_indices)
                best_solutions = population(valid_indices, :);
                best_fitness_values = fitness(valid_indices, :);
                for vi = 1:size(best_fitness_values,1)
                    if (best_fitness_values(vi,1)*w(1)+best_fitness_values(vi,2)*w(2)+best_fitness_values(vi,3)*w(3)+best_fitness_values(vi,4)*w(4)) ...
                            <= (global_best(1)*w(1)+global_best(2)*w(2)+global_best(3)*w(3)+global_best(4)*w(4))
                        global_solution = best_solutions(vi, :);
                        [global_best, best_Assoc, best_SINR, best_SINR_max, best_TP, best_DD] = objective_functions(global_solution);
                        fmin = global_best(1)*w(1)+global_best(2)*w(2)+global_best(3)*w(3)+global_best(4)*w(4);
                    end
                end
            end
        end

        gen = gen + 1;

        % Display progress
        %disp(["Generation: " num2str(gen) " Best Rank-1 Size: " num2str(length(fronts{1}))]);
    end

                %     best_Assoc = Assoc(valid_indices);
                % best_SINR = Mean_SINR(valid_indices);
                % best_DD = Mean_TP_eMBB(valid_indices);
                % best_SINR_max = SINR_max(:, valid_indices);

    final_metrics = [best_SINR, best_TP, best_DD, best_Assoc];

    Simulation_time_NSGA2 = toc;
    
    if plots == true

        X_CS = [];
        Y_CS = [];
        
        for i = 1:2:length(global_solution(1,:))
        
            X_CS = horzcat(X_CS,global_solution(1,i));
            Y_CS = horzcat(Y_CS,global_solution(1,i+1));
        
        end
        
        CEN_CS = transpose([X_CS; Y_CS]);
        
        for i = 1:1:length(global_solution(1,:))/2
            Rad(1,i) = R_limit;
        end
        
        figure
            hold on
            colorbar
            title('SINR (dB) by NSGA-II');
            scatter(Users(:,1),Users(:,2), 40, best_SINR_max(:,1), "filled", "MarkerFaceAlpha", 0.85)
            xlabel('x (km)') 
            ylabel('y (km)')
            viscircles(CEN_CS, Rad, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'red');
            scatter(X_CS(1,:),Y_CS(1,:), 400, 'pentagram', 'MarkerEdgeColor', 'red', 'linewidth', 1.5)
            rectangle('Position', [0 0 X Y])
            axis([-X/10 X+X/10 -Y/10 Y+Y/10])
            daspect([1 1 1])
            xlabel('x (m)') 
            ylabel('y (m)')
            hold off
    end

    function [z, Assoc, Mean_SINR, SINR_max, Mean_TP_eMBB, Mean_DD] = objective_functions(u)
    
    DronePop = zeros(2,N_UAV);

    for e = 1:2:2*N_UAV
        DronePop(:,(e+1)/2) = [u(1,e);u(1,e+1)];
    end
    
    N_users = length(Users);
    
    f = 3.5 * 10^9;
    velc = 299792458;
    ZetaLOS = 1;
    ZetaNLOS = 20;
    B = [50*10^6 20*10^6 10*10^6];
    q = -174 + 10*log10(B);
    
    Gt = 3;
    Gr = 0;

    for b = 1:1:N_UAV
        for k = 1:1:N_users
    
            D(k,b) = sqrt(((Users(k,1)-DronePop(1,b)).^2) + (Users(k,2)-DronePop(2,(b))).^2 + h.^2);
            R(k,b) = sqrt(abs(((Users(k,1)-DronePop(1,b)).^2) + (Users(k,2)-DronePop(2,(b))).^2));
            Demand_Density(k,b) = Demand(k)./(D(k,b).^2);
    
             switch Demand(k)
        
                 case DataRate(1)
                     User_Service(k) = 1;
                     
                 case DataRate(2)
                     User_Service(k) = 2;
        
                 case DataRate(3)
                     User_Service(k) = 3;
        
                 otherwise
                     error("Out of service");
             end
    
        end
    end

    theta = atan(h./R);

    values = [alpha, beta];

    Z = (alpha*exp(-beta*((180/pi).*theta - alpha)));

    PL = 20*log10((4*pi*f*D./velc))+((ZetaLOS+Z.*ZetaNLOS)./(1+Z));

    Pr_user = Pt - PL + Gt + Gr;

    Pr_Linear = (10.^((Pr_user-30)./10));

    % SINR Estimation

    for k = 1:1:N_UAV

        for j = 1:1:N_users

            inter = 0;

                for m = 1:1:N_UAV

                    if (k ~= m)
    
                        inter = inter + Pr_Linear(j,m);
    
                    end
                end

            Interference_lin(j,k) = inter;

            Interference(j,k) = 10*log10(inter);

            switch Demand(j)

                case DataRate(1)
                if R(j,k) > R_limit
                    SINR_lin(j,k) = 0;
                    TP(j,k) = 0;
                else
                    SINR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(1)-30)/10)) + inter);
                    SNR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(1)-30)/10)));
                    TP(j,k) = (10^-6) * (B(1))*log2(1+SNR_lin(j,k));
                end

                SINR(j,k) = 10.*log10(SINR_lin(j,k));

                case DataRate(2)

                if R(j,k) > R_limit
                    SINR_lin(j,k) = 0;
                    TP(j,k) = 0;
                else
                    SINR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(2)-30)/10)) + inter);
                    SNR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(2)-30)/10)));
                    TP(j,k) = (10^-6) * (B(2))*log2(1+SNR_lin(j,k));
                end

                SINR(j,k) = 10.*log10(SINR_lin(j,k));

                case DataRate(3)

                if R(j,k) > R_limit
                    SINR_lin(j,k) = 0;
                    TP(j,k) = 0;
                else
                    SINR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(3)-30)/10)) + inter);
                    SNR_lin(j,k) = Pr_Linear(j,k) / ((10^((q(3)-30)/10)));
                    TP(j,k) = (10^-6) * (B(3))*log2(1+SNR_lin(j,k));
                end

                SINR(j,k) = 10.*log10(SINR_lin(j,k));

            end

        end

    end

    Assoc = 0;
    TP_eMBB = 0;
    eMBB_count = 0;
    
    % for j = 1:1:size(SINR, 3)
    % 
    %     SINR_positioning(:,j) = SINR(:,j);
    %     TP_positioning(:,j) = TP(:,j);
    %     R_positioning(:,j) = R(:,j);
    % 
    % end
    
    Demand_max = max(Demand_Density,[],2);
    SINR_lin_max = max(SINR_lin,[],2);
    SINR_max = 10.*log10(SINR_lin_max);
    TP_max = max(TP,[],2);
    R_min = min(R,[],2);
    
    for k = 1:1:size(Users,1)
    
         if SINR_max(k,1) >= -10 && TP_max(k,1) > 0 && R_min(k,1) <= R_limit
            
             Assoc = Assoc + 1;
            
            if User_Service(k) == 1
                TP_eMBB = TP_eMBB + TP_max(k);
                 eMBB_count = eMBB_count + 1;
            end

        end
    
    end
    
    Mean_SINR_lin = mean(SINR_lin_max,1);
    
    Mean_SINR = 10*log10(Mean_SINR_lin);

    Mean_TP_eMBB = TP_eMBB/eMBB_count;
    
    Mean_DD = mean(Demand_max,1);
    
    % Cost Functions
    
    z1 = -Mean_SINR/normal_values(1);
    
    z2 = -Mean_TP_eMBB/normal_values(2);
    
    z3 = -Mean_DD/normal_values(3);
    
    z4 = -Assoc/normal_values(4);

    z = [z1, z2, z3, z4];
end

function offspring = crossover(parents, crossover_prob, lb, ub)
    num_parents = size(parents, 1);
    num_variables = size(parents, 2);
    offspring = parents;
    eta_c = 10;  % SBX parameter (higher = less diversity)

    for i = 1:2:num_parents-1
        if rand < crossover_prob
            p1 = parents(i, :);
            p2 = parents(i+1, :);
            gamma = zeros(1, num_variables);
            for j = 1:num_variables
                if rand < 0.5
                    gamma(j) = (2 * rand)^(1 / (eta_c + 1));
                else
                    gamma(j) = (1 / (2 * (1 - rand)))^(1 / (eta_c + 1));
                end
            end
            offspring(i, :) = 0.5 * ((1 + gamma) .* p1 + (1 - gamma) .* p2);
            offspring(i+1, :) = 0.5 * ((1 - gamma) .* p1 + (1 + gamma) .* p2);
        end
    end
end


function mutated_offspring = mutation(offspring, mutation_prob, lb, ub)
    num_offspring = size(offspring, 1);
    num_variables = size(offspring, 2);
    mutated_offspring = offspring;

    for i = 1:num_offspring
        for j = 1:num_variables
            if rand < mutation_prob
                mutated_offspring(i, j) = lb(j) + rand * (ub(j) - lb(j));
            end
        end
    end
end

function selected = tournament_selection(population, ranks, crowding_distances)
    n = size(population, 1);
    selected = zeros(size(population));

    for i = 1:n
        a = randi(n);
        b = randi(n);

        if (ranks(a) < ranks(b)) || (ranks(a) == ranks(b) && crowding_distances(a) > crowding_distances(b))
            selected(i, :) = population(a, :);
        else
            selected(i, :) = population(b, :);
        end
    end
end


function [fronts, ranks] = non_dominated_sorting(fitness)
    n = size(fitness, 1);
    ranks = zeros(n, 1);
    domination_count = zeros(n, 1);
    dominated_solutions = cell(n, 1);
    fronts = cell(1, n);
    
    for i = 1:n
        dominated_solutions{i} = [];
    end
    first_front = [];

    % Non-dominated sorting
    for i = 1:n
        for j = i+1:n  
            if all(fitness(i, :) <= fitness(j, :)) && any(fitness(i, :) < fitness(j, :))
                dominated_solutions{i} = [dominated_solutions{i}, j];
                domination_count(j) = domination_count(j) + 1;
            elseif all(fitness(j, :) <= fitness(i, :)) && any(fitness(j, :) < fitness(i, :))
                dominated_solutions{j} = [dominated_solutions{j}, i];
                domination_count(i) = domination_count(i) + 1;
            end
        end
        if domination_count(i) == 0
            ranks(i) = 1;
            first_front = [first_front, i];
        end
    end

    fronts{1} = first_front;
    front_idx = 1;

    while ~isempty(fronts{front_idx})
        next_front = [];
        for i = fronts{front_idx}
            for j = dominated_solutions{i}
                domination_count(j) = domination_count(j) - 1;
                if domination_count(j) == 0
                    ranks(j) = front_idx + 1;
                    next_front = [next_front, j];
                end
            end
        end
        if ~isempty(next_front)
            fronts{front_idx + 1} = next_front;
        end
        front_idx = front_idx + 1;
    end
end

function distances = crowding_distance(fitness, fronts)
    num_solutions = size(fitness, 1);
    distances = zeros(num_solutions, 1);

    for f = 1:length(fronts)
        front = fronts{f};
        if length(front) > 2
            for m = 1:size(fitness, 2)
                [~, sorted_idx] = sort(fitness(front, m));
                distances(front(sorted_idx(1))) = inf;
                distances(front(sorted_idx(end))) = inf;
                for i = 2:length(sorted_idx)-1
                    distances(front(sorted_idx(i))) = distances(front(sorted_idx(i))) + ...
                        (fitness(front(sorted_idx(i+1)), m) - fitness(front(sorted_idx(i-1)), m));
                end
            end
        end
    end
end

    function [new_population, new_fitness] = select_next_generation(combined_population, combined_fitness, n)
    % Perform non-dominated sorting
    [fronts, ranks] = non_dominated_sorting(combined_fitness);
    
    % Initialize new population
    new_population = [];
    new_fitness = [];
    
    front_idx = 1;
    
    % Add complete Pareto fronts until we reach the population limit
    while front_idx <= length(fronts) && ~isempty(fronts{front_idx}) ...
        && (size(new_population, 1) + length(fronts{front_idx}) <= n)

        new_population = [new_population; combined_population(fronts{front_idx}, :)];
        new_fitness = [new_fitness; combined_fitness(fronts{front_idx}, :)];
        front_idx = front_idx + 1;
    end
    
    % If the next front exceeds the population limit, use crowding distance
    remaining_spots = n - size(new_population, 1);
    if remaining_spots > 0
        % Compute crowding distances for the next front
        distances = crowding_distance(combined_fitness, fronts);
        
        % Get individuals from the front, sorted by crowding distance (descending)
        next_front = fronts{front_idx};
        [~, sorted_indices] = sort(distances(next_front), 'descend');
        selected_indices = sorted_indices(1:remaining_spots);
        
        % Add the best individuals based on crowding distance
        new_population = [new_population; combined_population(next_front(selected_indices), :)];
        new_fitness = [new_fitness; combined_fitness(next_front(selected_indices), :)];
    end
end

end