%% Network parameters
Nx = 8;
Ny = 8;
N  = Nx * Ny;       % 900 neurons

% Neuron type:
% 1 = PY (pyramidal)
% 2 = IN (interneuron)

type = ones(N,1);   % initially all PY
type(5:5:N) = 2;    % every fifth neuron is IN

% Coordinates
[x,y] = meshgrid(1:Nx,1:Ny);

x = x(:);
y = y(:);

%% Connection probabilities

p_PY_PY = 0.4;
p_PY_IN = 0.3;
p_IN_PY = 0.6;
p_IN_IN = 0.7;

LF_PY_PY = 10;
LF_PY_IN = 10;
LF_IN_PY = 20;

% Connectivity matrix
% C(i,j) = 1 means neuron j -> neuron i
C = zeros(N,N);

%% Generate connections

for pre = 1:N

    for post = 1:N

        % Don't connect a neuron to itself
        if pre == post
            continue
        end

        % Spatial distance
        dx = abs(x(pre) - x(post));
        dy = abs(y(pre) - y(post));

        % Determine connection type
        preType  = type(pre);
        postType = type(post);

        % PY -> PY
        if preType == 1 && postType == 1

            LF = LF_PY_PY;
            p  = p_PY_PY;

        % PY -> IN
        elseif preType == 1 && postType == 2

            LF = LF_PY_IN;
            p  = p_PY_IN;

        % IN -> PY
        elseif preType == 2 && postType == 1

            LF = LF_IN_PY;
            p  = p_IN_PY;

        % IN -> IN
        elseif preType == 2 && postType == 2

            % Paper gives probability but no separate LF
            LF = 10;
            p  = p_IN_IN;
        end

        % Is post neuron inside the square footprint?
        if dx <= LF && dy <= LF

            % Randomly establish connection
            if rand < p
                C(post,pre) = 1;
            end

        end
    end
end
fprintf('Total neurons = %d\n',N);
fprintf('PY neurons    = %d\n',sum(type==1));
fprintf('IN neurons    = %d\n',sum(type==2));
C