function [grid]=makegrid(mpar,grid)


%% Quadruble Log Grid
m_min = 0; %natural borrowing limit (adjust according to eq interest rate)
m_max = 150; % 150
grid.m = exp(exp(linspace(0,log(log(m_max - m_min+1)+1),mpar.nm))-1)-1+m_min;

% grid.m(abs(grid.m)==min(abs(grid.m)))=0;
