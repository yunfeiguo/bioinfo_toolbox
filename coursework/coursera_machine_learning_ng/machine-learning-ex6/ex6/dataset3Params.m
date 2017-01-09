function [C, sigma] = dataset3Params(X, y, Xval, yval)
%DATASET3PARAMS returns your choice of C and sigma for Part 3 of the exercise
%where you select the optimal (C, sigma) learning parameters to use for SVM
%with RBF kernel
%   [C, sigma] = DATASET3PARAMS(X, y, Xval, yval) returns your choice of C and 
%   sigma. You should complete this function to return the optimal C and 
%   sigma based on a cross-validation set.
%

% You need to return the following variables correctly.
C = 30;
sigma = 30;

% ====================== YOUR CODE HERE ======================
% Instructions: Fill in this function to return the optimal C and sigma
%               learning parameters found using the cross validation set.
%               You can use svmPredict to predict the labels on the cross
%               validation set. For example, 
%                   predictions = svmPredict(model, Xval);
%               will return the predictions on the cross validation set.
%
%  Note: You can compute the prediction error using 
%        mean(double(predictions ~= yval))
%
C_min = 1;
sigma_min = 1;
min_error = 1;
if true
    for C = [0.01,0.03,0.1,0.3,1,3,10,30]
	for sigma = [0.01,0.03,0.1,0.3,1,3,10,30]
	    % Train the SVM
	    model= svmTrain(X, y, C, @(x1, x2) gaussianKernel(x1, x2, sigma));
	    current_error = mean(double(yval ~= svmPredict(model, Xval)));
	    if current_error < min_error
		min_error = current_error;
		C_min = C;
		sigma_min = sigma;
	    end

	fprintf('C %f and sigma %f: %f\n', C, sigma, current_error);
end
end
end
fprintf('C min %f and sigma min %f: %f\n', C_min, sigma_min, min_error);

C = C_min;
sigma = sigma_min;




% =========================================================================

end
