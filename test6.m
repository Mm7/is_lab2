
% Pick some random epsilon and delta for each of them do the following
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

    % Collect some statistics:
    %   `p_leg`: probability of mis-decode by the legitimate receiver.
    p_leg = 0;
    cnt = 10000;

    for i = 1:cnt
        % Generate a random message.
        u = randi(8) - 1;

        % Encode, send through BSC and decode.
        x = uniform_encode(u);
        [y, ~] = bsc_channel(x, epsilon, 0);
        u_y = uniform_decode(y);

        % Save the results of this simulation.
        if u ~= u_y
            p_leg = p_leg + 1;
        end
    end

    % Normalize the counts to get empirical probabilities.
    p_leg = p_leg / cnt;

    % Observation: perfect reliability is not guaranteed anymore.
    fprintf("\t\tProbabiliy of mis-decoding the message by the legitimate receiver: %.1f%%\n", p_leg*100);

    %% Reproduce the simulation of task4 with a BSC.

    fprintf("\t2) Reproducing task4 with a BSC\n");

    % p_z|u
    p_z_u = zeros(8, 128);

    % # of samples to generate to estimate a probability value.
    cnt = 10000;

    for u = 0:7   
        % Simulation...
        for i = 1:cnt
            z = eavesdropper(u, delta);
            p_z_u(u+1, z+1) = p_z_u(u+1, z+1) + 1;
        end
    end

    p_z_u = p_z_u / cnt;

    % Plot p_z|u.
    figure();
    for i = 1:8
        scatter(1:128, p_z_u(i,:) * 100);
        hold on;
    end

    title(sprintf("p_{z|u} - epsilon: %.2f delta: %.2f", epsilon, delta));
    subtitle("Each line is a different message u");
    grid on;
    xlim([0, 128]);
    ylabel("Probability [%]");
    xlabel("Codeword [z]");
    legend(["0b000","0b001","0b010","0b011","0b100","0b101","0b110","0b111"]);

    % Compute some PMDs and I(u, z).
    p_u = ones(8, 1) * (1/8);
    p_uz = p_z_u .* (ones(8, 128) * (1/8));
    p_z = sum(p_uz, 1)';

    I_uz = 0;

    for u = 1:8
        for z = 1:128
            if p_uz(u, z) ~= 0
                I_uz = I_uz + p_uz(u, z) * log2(p_uz(u, z) / (p_u(u) * p_z(z)));
            end
        end
    end

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
    
    % Unfortunately `p_leg` computed before is not the same as `p_max`. So
    % recompute it properly.
    p_max = 0;
    cnt = 10000;

    for u = 0:6
        % P(u = u_y | u) 
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
    
    % The mutual information on the other hand should be correct.
    half_sq_I = 0.5 * sqrt(I_uz);
    
    % Compute the `e` security parameter. If the upper bound is higher than
    % 1 then throw it and use simply 1 (which is the trivial upper bound).
    e = min(p_max + half_sq_I, 1);
    
    fprintf("\t\tMax error probability: %.2f\n", p_max);
    fprintf("\t\tHalf sqrt of I: %.2f\n", half_sq_I);
    fprintf("\t\tSecurity parameter upper bound: %.2f\n", e);
end

%% Helper functions

% Compute the entropy of a random variable given its PMD.
function entropy = entropy(v)
    lg = -log2(v);
    lg(isinf(lg)) = 0;
    entropy = sum(lg .* v, 'all');
end

% Chain the encoding + wiretrap channel.
function z = eavesdropper(u, delta)
    x = uniform_encode(u);
    [~,z] = bsc_channel(x, 0, delta);
end
