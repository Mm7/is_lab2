%% Plot the PMD p_z|u

% p_z|u
p_z_u = zeros(8, 128);

% # of samples to generate to estimate a probability value.
COUNT = 10000;

for u = 0:7    
    % Simulation...
    for i = 1:COUNT
        z = eavesdropper(u);
        
        p_z_u(u+1, z+1) = p_z_u(u+1, z+1) + 1;
    end
end

p_z_u = p_z_u / COUNT;

% Plot p_z|u.
for i = 1:8
    plot(1:128, p_z_u(i,:) * 100);
    hold on;
end

title("p_{z|u}");
subtitle("Each line is a different message u");
grid on;
xlim([0, 128]);
ylim([0, 1.8]);
ylabel("Probability [%]");
xlabel("Word [z]");
legend(["0b000","0b001","0b010","0b011","0b100","0b101","0b110","0b111"]);

%% Compute some PMDs and I(u, z).
p_u = ones(8, 1) * (1/8);
p_uz = p_z_u .* (ones(8, 128) * (1/8));
p_z = sum(p_uz, 1)';

I_uz = 0;

for u = 1:8
    for z = 1:128
        I_uz = I_uz + p_uz(u, z) * log2(p_uz(u, z) / (p_u(u) * p_z(z)));
    end
end

fprintf("Entropy of u = %f\n", entropy(p_u));
fprintf("Entropy of z = %f\n", entropy(p_z));
fprintf("Joint entropy of u,z = %f\n", entropy(p_uz));
fprintf("Mutual information I(u,z) = %.3f\n", I_uz);

%% Helper functions

% Compute the entropy of a random variable given its PMD.
function entropy = entropy(v)
    entropy = sum(-log2(v) .* v, 'all');
end

% Chain the encoding + wiretrap channel.
function z = eavesdropper(u)
    x = uniform_encode(u);
    
    % Convert `x` to a bitarray.
    xi = zeros(1,7);
    for j = 1:7
        if bitget(x, 8-j) > 0
            xi(j) = 1;
        end
    end
    
    [~,zi] = wiretrap_channel(xi);

    % Convert `z` to an integer.
    z = 0;
    for j = 0:6
        if zi(1,7-j) > 0
            z = bitset(z, j+1);
        end
    end
end
