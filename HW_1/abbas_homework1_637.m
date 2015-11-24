% This is a demonstration program to get you started in Matlab: To Create a
% Program or 'Script': click on the 'New' icon on the top left side of the
% Matlab home page. Click on the 'Save icon and name the program
% 'first_program.m' Back in the Script type the following code:
%% Line Plots
%%
% To create two-dimensional line plots, use the |plot| function. For
% example, plot the value of the sine function from 0 to $2\pi$:
x = 0:pi/100:2*pi; %Create a. vector 'x' with values %between 0 and 2*pi * in increments of pi/100.
y = sin(x); % Create a vector 'y' using the values of the x vector as * arguments of a sine function.
plot(x,y) % plot a 2-D graph with y as the domain and x. the range. You can label the axes and add a title.
xlabel('x')
ylabel('sin(x)')
title('Plot of the Sine Function')

% By adding a third input argument to the |plot| function, you can plot the
% % same variables using a dashed red, indicated by 'r', line.
plot(x,y,'c--')
%%
% The |'r--'| string is a _line specification_. Each specification can
% include characters for the line color, style, and marker. A style can be
% a dashed '--' or solid '-' line. A marker is a symbol that appears at
% each plotted data point, such as a, |+|, |0|, or |*|. For example,
% |’g:*’| requests a dotted green line with |*| markers. Notice that the
% titles and labels that you defined for the first plot are no longer in
% the current figure_ window. By default, MATLAB(R) clears the figure each
% time you call a plotting function, resetting the axes and other elements
% to prepare the new plot.
%%
% To add plots to an existing figure, use |hold|.
x = 0:pi/100:2*pi;
y = sin(x);
plot (x, y)
hold on
y2 = cos(x);
plot(x,y2,':')
legend('sin','cos')
% Until you use |hold off| or close the window, all plots appear in the
% current figure window.
%%
% Copyright 2014, The MathWorks Inc.