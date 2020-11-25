
% Run in fast mode to generate all the plots quickly at the expense of a
% lower granularity and accuracy.
global fast;
fast = true;

%% Run some simulations for different pairs of epsilon and delta.

% Pick some random epsilon and delta, for each of them do the following
% steps..
for k = 1:4
    % Random parameters..
    epsilon = rand() / 2;
    delta = rand() / 2;

    fprintf(">>> BSC parameters: epsilon: %.2f, delta: %.2f\n", epsilon, delta);

    %% Reproduce the simulation of task3 with a BSC.

    fprintf("\t1) Reproducing task3 with a BSC\n");

    % In this case we expect to see some mis-decoding for the legitimate
    % receiver.
    p_leg = decoding_err_prob(epsilon);

    % Observation: perfect reliability is not guaranteed anymore.
    fprintf("\t\tProbabiliy of mis-decoding the message by the legitimate receiver: %.1f%%\n", p_leg*100);

    %% Reproduce the simulation of task4 with a BSC.

    fprintf("\t2) Reproducing task4 with a BSC\n");
    
    % Compute all the required values.
    [p_z_u, p_uz, p_u, p_z, I_uz] = compute_probs(delta);

    % Plot p_z|u.
    figure();
    for i = 1:8
        scatter(1:128, p_z_u(i,:) * 100);
        hold on;
    end

    title(sprintf("p_{z|u} - epsilon: %.2f delta: %.2f", epsilon, delta));
    subtitle("Each color corresponds to a different message u");
    grid on;
    xlim([0, 128]);
    ylabel("Probability [%]");
    xlabel("Word [z]");
    legend(["0b000","0b001","0b010","0b011","0b100","0b101","0b110","0b111"]);

    % Print some values.
    fprintf("\t\tEntropy of u = %f\n", entropy(p_u));
    fprintf("\t\tEntropy of z = %f\n", entropy(p_z));
    fprintf("\t\tJoint entropy of u,z = %f\n", entropy(p_uz));
    fprintf("\t\tMutual information I(u,z) = %.3f\n", I_uz);
    
    %% Compute the e-security of the mechanism.
    
    % Let:
    %   p_max = maximum (over d in M) of P(u != u_recovered | u = d )
    %           conditional decoding error probability.
    %   I(u,z) = mutual information between `u` and `z`.
    % An upper bound on the e-security of the mechanism is:
    %   p_e + sqrt(I(u,z))/2.
    
    fprintf("\t3) Computing the e-security of the mechanism\n");
    
    e = e_security(epsilon, I_uz);

    fprintf("\t\tSecurity parameter upper bound: %.2f\n", e);
end

%% Plot 2.3.4: decoding error probability versus epsilon.

fprintf(">>> Plot 2.3.4: decoding error probability vs epsilon\n");

granularity = 0.001 * max(fast * 10, 1);
epsilon_vals = 0:granularity:0.5;
prob_err = zeros(length(epsilon_vals), 1);

for i = 1:length(epsilon_vals)
    prob_err(i) = decoding_err_prob(epsilon_vals(i));
end

figure();

plot(epsilon_vals, prob_err*100);

title("Decoding error probabiliy P(u != u_{decoded}) vs espilon");
xlabel("Epsilon");
ylabel("Decoding error probability [%]");
grid on;

%% Plot 2.3.5: Mutual information versus delta.

fprintf(">>> Plot 2.3.5: mutual information vs I(u,z)\n");

granularity = 0.001 * max(fast * 10, 1);
delta_vals = 0:granularity:0.5;
mut_inf = zeros(length(delta_vals), 1);

for i = 1:length(delta_vals)
    [~, ~, ~, ~, I_uz] = compute_probs(delta_vals(i));
    mut_inf(i) = I_uz;
end

figure();

plot(delta_vals, mut_inf);

title("Mutual information I(u,z) vs delta");
xlabel("Delta");
ylabel("I(u,z)");
grid on;

%% Plot 2.3.6: Upper bound on security of the mechanism versus (epsilon, delta).

fprintf(">>> Plot 2.3.6: upper bound on the security of the mechanism vs (epsilon, delta)\n");

granularity = 0.01 * max(fast * 10, 1);
epsilon_vals = 0:granularity:0.5;
delta_vals = 0:granularity:0.5;

e_sec = zeros(length(epsilon_vals), length(delta_vals));

for i = 1:length(epsilon_vals)
    for j = 1:length(delta_vals)
        [~, ~, ~, ~, I_uz] = compute_probs(delta_vals(j));
        e_sec(i, j) = e_security(epsilon_vals(i), I_uz);
    end
end

figure();

[C,h] = contour(epsilon_vals, delta_vals, e_sec);
clabel(C, h);

grid on;
xlabel("Delta");
ylabel("Epsilon");

title("Upper bound on security parameter vs (epsilon, delta)");

%% Helper functions

% Compute the decoding error probability over a BSC, given the epsilon
% parameter.
function p = decoding_err_prob(epsilon)
    global fast;
    p = 0;
    cnt = 1000 * max(~fast * 10, 1);

    for i = 1:cnt
        % Generate a random message.
        u = randi(8) - 1;

        % Encode, send through BSC and decode.
        x = uniform_encode(u);
        [y, ~] = bsc_channel(x, epsilon, 0);
        u_y = uniform_decode(y);

        % Save the results of this iteration.
        if u ~= u_y
            p = p + 1;
        end
    end

    % Normalize the counts to get empirical probabilities.
    p = p / cnt;
end

% Given the parameters epsilon delta of a BSC, run some simulations and
% compute the following empirical probabilities:
%   p_z_u: P(u|z)
%   p_uz:  P(u,z)
%   p_u:   P(u) - should be an input value. For now let's hardcode 1/8 for
%          all the outcomes.
%   p_z:   P(z)
%   I_uz:  I(u,z)
function [p_z_u, p_uz, p_u, p_z, I_uz] = compute_probs(delta)
    global fast;
    
    % Gather the statistics on the joint distribution of u and z.
    p_uz = zeros(8, 128);
    cnt = 1000 * max(~fast * 10, 1);

    for i = 1:cnt
        u = randi(8)-1;
        z = eavesdropper(u, delta);
        
        p_uz(u+1, z+1) = p_uz(u+1, z+1) + 1;
    end
    
    p_uz = p_uz / cnt;

    % Compute some PMDs and I(u, z).    
    p_u = sum(p_uz, 2);
    p_z = sum(p_uz, 1)';
    p_z_u = p_uz ./ (p_u * ones(1,128));

    I_uz = 0;

    for u = 1:8
        for z = 1:128
            if p_uz(u, z) ~= 0
                I_uz = I_uz + p_uz(u, z) * log2(p_uz(u, z) / (p_u(u) * p_z(z)));
            end
        end
    end
end

% Compute an upper bound on the security parameter of the mechanism given
% the `epsilon` parameter of the BSC and the mutual information I(u,z).
function e = e_security(epsilon, I_uz)
    global fast;

    % Let:
    %   p_max = maximum (over d in M) of P(u != u_recovered | u = d )
    %           conditional decoding error probability.
    %   I(u,z) = mutual information between `u` and `z`.
    % An upper bound on the e-security of the mechanism is:
    %   p_max + sqrt(I(u,z))/2.

    % First, compute `p_max` by the means of a simulation.
    p_max = 0;
    cnt = 1000 * max(~fast * 10, 1);

    for u = 0:6
        % Letting `u` the transmitted message and `u_y` the decoded message
        % (called `u_recovered` above) by the legitimate receiver, `p`
        % collects the statistics of:
        %   P(u != u_y | u)
        p = 0;
        
        for i = 1:cnt
            % Encode, send through BSC and decode.
            x = uniform_encode(u);
            [y, ~] = bsc_channel(x, epsilon, 0);
            u_y = uniform_decode(y);

            % Save the results of this simulation.
            if u ~= u_y
                p = p + 1;
            end
        end
    
        % Normalize the counts to get empirical probabilities.
        p = p / cnt;
        
        % Take the maximum.
        p_max = max([p, p_max]);
    end

    half_sq_I = 0.5 * sqrt(I_uz);
    
    % Compute the `e` security parameter. If the upper bound is higher than
    % 1 then throw it and use simply 1 (which is the trivial upper bound).
    e = min(p_max + half_sq_I, 1);
end

% Compute the entropy of a random variable given its PMD.
function entropy = entropy(v)
    lg = -log2(v);
    lg(isinf(lg)) = 0;
    entropy = sum(lg .* v, 'all');
end

% Chain the encoding + BSC channel.
function z = eavesdropper(u, delta)
    x = uniform_encode(u);
    [~,z] = bsc_channel(x, 0, delta);
end
