function [] = draw_fixation(p,w)
% parameters for drawing fixation circle given pre-set params
% re-used many many times in WM tasks w circular aperture
Screen('DrawDots',w,[0;0],p.fix_size_mult*p.fix_size_out*p.ppd*p.fix_size_constant+p.fix_pen,p.dim_amt*p.fix_color,p.center,2);
Screen('DrawDots',w,[0;0],p.fix_size_mult*p.fix_size_out*p.ppd*p.fix_size_constant-p.fix_pen,p.bg_color,p.center,2);
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*p.fix_size_constant,p.dim_amt*p.fix_color,p.center,2);
end

