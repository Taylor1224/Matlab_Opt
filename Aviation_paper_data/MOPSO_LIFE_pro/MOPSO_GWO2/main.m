%% main 

% Parameters
num_particles = 50; % Number of particles
max_iter = 50; % Maximum number of iterations
dim = 4; % Dimension of the problem
lb = [996, 646, 225, 490]; % Lower bound
ub = [1036, 675, 235, 510]; % Upper bound

% PSO Parameters
w = 0.5; % Inertia weight
c1 = 1.5; % Cognitive coefficient
c2 = 1.5; % Social coefficient

% Initialize particles
particles = lb + (ub - lb) .* rand(num_particles, dim); % Positions
velocities = zeros(num_particles, dim); % Velocities
personal_best_positions = particles; % Personal best positions
personal_best_values = inf(num_particles, 2); % Personal best values (for multi-objective)

% Initialize alpha, beta, delta wolves (best, second best, third best)
alpha_position = zeros(1, dim);
beta_position = zeros(1, dim);
delta_position = zeros(1, dim);
alpha_score = inf;
beta_score = inf;
delta_score = inf;

% Initialize archive for Pareto front
archive_positions = [];
archive_values = [];

% Main loop
for iter = 1:max_iter
    for i = 1:num_particles
        % Fitness evaluation (multi-objective)
        fitness = objective_function(particles(i, :));
        
        % Update personal best
        if dominates(fitness, personal_best_values(i, :))
            personal_best_values(i, :) = fitness;
            personal_best_positions(i, :) = particles(i, :);
        end
        
        % Update alpha, beta, delta
        if dominates(fitness, alpha_score)
            delta_score = beta_score;
            delta_position = beta_position;
            
            beta_score = alpha_score;
            beta_position = alpha_position;
            
            alpha_score = fitness;
            alpha_position = particles(i, :);
        elseif dominates(fitness, beta_score)
            delta_score = beta_score;
            delta_position = beta_position;
            
            beta_score = fitness;
            beta_position = particles(i, :);
        elseif dominates(fitness, delta_score)
            delta_score = fitness;
            delta_position = particles(i, :);
        end
    end

    % Update the position of particles using GWO + PSO
    for i = 1:num_particles
        % GWO update
        A1 = 2 * rand(dim, 1) - 1;
        C1 = 2 * rand(dim, 1);
        D_alpha = abs(C1 .* alpha_position' - particles(i, :)');
        X1 = alpha_position' - A1 .* D_alpha;
        
        A2 = 2 * rand(dim, 1) - 1;
        C2 = 2 * rand(dim, 1);
        D_beta = abs(C2 .* beta_position' - particles(i, :)');
        X2 = beta_position' - A2 .* D_beta;
        
        A3 = 2 * rand(dim, 1) - 1;
        C3 = 2 * rand(dim, 1);
        D_delta = abs(C3 .* delta_position' - particles(i, :)');
        X3 = delta_position' - A3 .* D_delta;
        
        % New position based on GWO
        new_position = (X1 + X2 + X3) / 3;

        % PSO velocity update
        velocities(i, :) = w * velocities(i, :) + ...
                           c1 * rand * (personal_best_positions(i, :) - particles(i, :)) + ...
                           c2 * rand * (alpha_position - particles(i, :));
        % Position update
        particles(i, :) = new_position' + velocities(i, :);
        
        % Boundary control
        particles(i, :) = max(particles(i, :), lb);
        particles(i, :) = min(particles(i, :), ub);
    end

    % Update the Pareto archive
    for i = 1:num_particles
        fitness = objective_function(particles(i, :));
        archive_positions = [archive_positions; particles(i, :)];
        archive_values = [archive_values; fitness];
    end
    % Perform non-dominated sorting and update the archive
    [archive_positions, archive_values] = update_archive(archive_positions, archive_values);

    % Display iteration information
    disp(['Iteration ' num2str(iter) ': Best alpha fitness = ' num2str(alpha_score)]);
end

% Output the Pareto front
disp('Pareto front solutions:');
disp(archive_positions);
disp('Objective values:');
disp(archive_values);




