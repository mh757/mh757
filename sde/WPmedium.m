classdef WPmedium < WP
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = protected)
        Xstart
        Xend
        DEBUG % verbosity for debug
    end

    methods
        function obj = WPmedium(N)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@WP(0, N, 1);
            obj.Xstart = 0;
            obj.Xend = 1;
            
        end

        function Wx = eval(obj, x)
            if (x ~= 0) && (x - floor(x) == 0)
                Wx = x * obj.path(end);
                idx = 1 + obj.L;
                dW_idx = obj.dW(idx-1);
            else
                xModT = x - floor(x); % on Torus
                idx = 1 + floor(xModT * obj.L);

                % Extrapolate
                d_xModT = xModT - (idx - 1) * obj.dt;
                dW_idx = obj.dW(idx); 
                slope = dW_idx / obj.dt;
                correction = floor(x) * obj.path(end);
                Wx = obj.path(idx) + slope * d_xModT + correction;
            end

            if (obj.DEBUG)
                fprintf('x = %f, idx = %d, dW(idx) = %f, W(x) = %f\n', ...
                    x, idx, dW_idx, Wx);
            end
        end

        function Wprime = d_dx(obj, x, epsilon)
            % d_dx - Take the discrete derivative at x.
            %
            Wprime = (obj.eval(x + epsilon) - obj.eval(x - epsilon)) / ...
                    (2 * epsilon);
        end

        function newpath = extend(obj)
            %extend - LPSX style extension of W
            %   Extends the medium

            N = obj.L;

            % original path will be of length N + 1
            newpath(N+1:2*N+1) = obj.path;
            
            % Extend to [-1, 0) following LPSX
            % i = 1 implies x = -1.
            for i = 1:N
                newpath(i) = newpath(N+i) - newpath(end);
            end
        end
    end
end