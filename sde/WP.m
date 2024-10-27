classdef WP
    % Wiener Process
    % ----------------------------------------------
    % Created using methods from the book
    % An Introduction to the Numerical Simulation of
    % Stochastic Differential Equations 
    % by Higham and Kloeden.
    % ----------------------------------------------
    % Author: Marcel Hudiani
    % 10/27/2024
    % ----------------------------------------------

    properties (Access = protected)
        L % Total number of sampling points from dt to T (length of array).
        T % Total simulation time.
        path % The realization of this process (array).
        dt % Discrete increments.
        dW % Discrete dW as a function of dt (array).
    end

    methods
        function obj = WP(L,T)
            % Construct an instance of this class
            obj.L = L;
            obj.T = T;
            obj.dt = T/L;

            % Brownian increments
            % W(t_{i+1}) = W(t_{i}) + dW
            obj.dW = sqrt(obj.dt) * randn(1, obj.L);
            obj.path = cumsum(obj.dW);

            % W(0) = 0 for standard BM
            obj.path = [0, obj.path];
        end

        function W = get_path(obj)
            W = obj.path;
        end

        function sampledTimes = get_sampling_times(obj)
            sampledTimes = 0:obj.dt:obj.T;
        end

        function [W, err] = integrate(obj, integrand, varargin)
            % Programming Exercise 4.2:
            % -------------------------------------------------------
            % Compute the Ito integral w.r.t this object where
            % the integrand is a function of time.
            % -------------------------------------------------------
            % TODO: Add Stratonovich option
            %--------------------------------------------------------

            % TODO: Assign optional arguments
            if ~isempty(varargin)
                if length(varargin) >= 1
                    isIto = false;
                    % compute Stratonovich correction
                end
            else
                isIto = true;
            end
            
            % Ito integral
            W = sum(integrand .* obj.dW);
            err = abs(W - 0.5 * (obj.path(end)^2 - obj.T));

            % TODO: Stratonovich correction
            if ~isIto
            end
        end

        function obj = fill(obj, varargin)
            % Programming Exercise 3.1:
            % -----------------------------------------------------------
            % How can we assign W(t_{i + alpha}) where alpha is in (0, 1)
            % to an existing W(t)?
            % W(t_{i + 1/2}) = Average{W(t_{i + 1}), W(t_i)} + N(0, dt/4)
            % ------------------------------------------------------------
            % TODO: Generalize for alpha != 0.5.
            % However, we can keep increasing the sampling frequency
            % by repeated calls of this method.
            % ------------------------------------------------------------
            
            % Initialization
            Wfill = zeros(1, obj.L);
            alpha = 0.5;

            % Assign optional arguments
            if ~isempty(varargin)
                if length(varargin) >= 1
                    alpha = varargin{1}; % First optional argument
                end
            end

            % Error handling: Alpha must be inside the interval (0, 1).
            if (alpha <= 0 || alpha >= 1)
                alpha = 0.5;
            end

            % Figure out the filling separately for code readability.
            for t = 1:obj.L
                Wfill(t) = 0.5*(obj.path(t) + obj.path(t+1)) + ...
                    0.5 * sqrt(obj.dt) * randn(1,1);
            end

            % Insert the required fill-in at the right time index
            for t = 1:obj.L
                insertAt = 2 * (t - 1) + 1;
                obj.path = [obj.path(1 : insertAt), ...
                    Wfill(t), ...
                    obj.path(insertAt + 1 : end)];
            end

            % Update this object's parameters (length of array)
            obj.L = length(obj.path) - 1;
            obj.dt = obj.dt * alpha;
            obj.dW = zeros(1,obj.L);
            for t = 1:obj.L
                obj.dW(t) = obj.path(t+1) - obj.path(t);
            end
        end
        % methods end
    end
    % class end
end