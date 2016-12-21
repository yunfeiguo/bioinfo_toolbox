function [J grad] = nnCostFunction(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
%NNCOSTFUNCTION Implements the neural network cost function for a two layer
%neural network which performs classification
%   [J grad] = NNCOSTFUNCTON(nn_params, hidden_layer_size, num_labels, ...
%   X, y, lambda) computes the cost and gradient of the neural network. The
%   parameters for the neural network are "unrolled" into the vector
%   nn_params and need to be converted back into the weight matrices. 
% 
%   The returned parameter grad should be a "unrolled" vector of the
%   partial derivatives of the neural network.
%

% Reshape nn_params back into the parameters Theta1 and Theta2, the weight matrices
% for our 2 layer neural network
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));

% Setup some useful variables
m = size(X, 1);
         
% You need to return the following variables correctly 
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));

% ====================== YOUR CODE HERE ======================
% Instructions: You should complete the code by working through the
%               following parts.
%
% Part 1: Feedforward the neural network and return the cost in the
%         variable J. After implementing Part 1, you can verify that your
%         cost function computation is correct by verifying the cost
%         computed in ex4.m
%
% Part 2: Implement the backpropagation algorithm to compute the gradients
%         Theta1_grad and Theta2_grad. You should return the partial derivatives of
%         the cost function with respect to Theta1 and Theta2 in Theta1_grad and
%         Theta2_grad, respectively. After implementing Part 2, you can check
%         that your implementation is correct by running checkNNGradients
%
%         Note: The vector y passed into the function is a vector of labels
%               containing values from 1..K. You need to map this vector into a 
%               binary vector of 1's and 0's to be used with the neural network
%               cost function.
%
%         Hint: We recommend implementing backpropagation using a for-loop
%               over the training examples if you are implementing it for the 
%               first time.
%
% Part 3: Implement regularization with the cost function and gradients.
%
%         Hint: You can implement this around the code for
%               backpropagation. That is, you can compute the gradients for
%               the regularization separately and then add them to Theta1_grad
%               and Theta2_grad from Part 2.
%

X = [ones(m,1),X];
y_binary = zeros(m,num_labels);
for i=1:num_labels
    y_binary(:,i) = y==i;
end
y = y_binary;

%function [J grad] = nnCostFunction(nn_params, ...
%                                   input_layer_size, ...
%                                   hidden_layer_size, ...
%                                   num_labels, ...
%                                   X, y, lambda)
%forward propagation
%X:m*(input_layer_size+1)
%Theta1:hidden_layer_size*(input_layer_size+1)
%Theta2:num_labels*(hidden_layer_size+1)
%y_binary:m*num_labels

%Theta1*X':hidden_layer_size*m
hidden_layer = sigmoid(Theta1*X')';
hidden_layer = [ones(m,1), hidden_layer];
output_layer = sigmoid(Theta2*hidden_layer')';

%why can we unroll? because all terms are independent of another
Theta1_no_bias = Theta1(:,2:end);
Theta2_no_bias = Theta2(:,2:end);
all_param_no_bias = [Theta1_no_bias(:);Theta2_no_bias(:)];
J = -1*(y(:)'*log(output_layer(:)) + (1-y(:))'*log(1-output_layer(:)))/m + lambda*(all_param_no_bias'*all_param_no_bias)/m/2;

%backward propagation
%for t = 1:m
%    delta_3 = output_layer(t,:)' - y(t,:)';
%    delta_2 = Theta2'*delta_3.*hidden_layer(t,:)'.*(1-hidden_layer(t,:)');
%    delta_2 = delta_2(2:end,:); %remove error for bias unit
%    Theta2_grad(:,2:end) += delta_3 * hidden_layer(t,2:end);
%    Theta2_grad(:,1) += delta_3 * hidden_layer(t,1);
%    Theta1_grad(:,2:end) += delta_2 * X(t,2:end);
%    Theta1_grad(:,1) += delta_2 * X(t,1);
%end
%same thing, but vectorized
delta_3 = output_layer - y; %m*k
delta_2 = (Theta2'*delta_3')'.*hidden_layer().*(1-hidden_layer);
delta_2 = delta_2(:,2:end); %remove error for bias unit


Theta2_grad = delta_3' * hidden_layer;
Theta1_grad = delta_2' * X;
Theta2_grad(:,2:end) += lambda*Theta2(:,2:end);
Theta1_grad(:,2:end) += lambda*Theta1(:,2:end);
Theta2_grad /= m;
Theta1_grad /= m;
% -------------------------------------------------------------

% =========================================================================

% Unroll gradients
grad = [Theta1_grad(:) ; Theta2_grad(:)];


end
