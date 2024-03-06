function [delta] = run_staircase(history,target_acc,lastdelta)
%run_staircase Based on history and target accuracy, run staircase
%   history is the last n trials' amount of correct trials
%   target_acc is the accuracy level we're aiming for (in percentage)
%   the output, delta, is the distance between target and test that should
%   be run next
%   Try a weighted up/down method where target_acc =
%   step_size_up/(step_size_up+step_size_down);
%   a.k.a. step_size_up = target_acc, step_size_down = 1-target_acc;
mindelta = 1;
maxdelta = 179;

ntrials = numel(history); 
movement_scalar = max([10*(ntrials<=20) 5]);% move 10 steps at a time until 20 trials of each condition have happened
% then move 5 steps at a time for the rest of the experiment
if history(end)==0
    step_size = target_acc; %make delta larger
else % make delta smaller, correct answer
    step_size = target_acc-1; %negative step size
end
step_size = step_size*movement_scalar; %run regular up/down for 20 trials of each condition
delta = lastdelta+step_size;
delta(delta<mindelta) = mindelta; delta(delta>maxdelta)=maxdelta;
end

