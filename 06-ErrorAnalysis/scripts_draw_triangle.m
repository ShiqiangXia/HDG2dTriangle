
start_x = 1.925;
start_y = -3;

end_x = 2;

slope = -7;

delta = abs(start_x-end_x);
end_y = start_y + slope * delta;

hold on;

% draw horizontal line
plot([start_x,end_x ],[start_y,start_y],...
    '-k','LineWidth',1,'HandleVisibility','off')

% draw vertical line
hold on;

plot([end_x,end_x],[start_y,end_y],...
    '-k','LineWidth',1,'HandleVisibility','off')

% draw line 3
plot([start_x, end_x ],[start_y,end_y],...
    '-k','LineWidth',1,'HandleVisibility','off')


