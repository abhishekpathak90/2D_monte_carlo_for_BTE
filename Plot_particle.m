movieLength = 5; % duration of movie in seconds
N =4000; % number of trajectories
for ii=1:N
    filename = ['Data_processing/' num2str(ii) '.txt'];
    Traj = load(filename);
    
    figure();
    hold on
    %removing zero rows
    %Traj(~any(Traj,2),:) = [];
    [row, col] = size(Traj);
    for i=1:row
        plot(Traj(i,2),Traj(i,3), 'ro');
        xlim([0 1e-7]);
        ylim([0 1e-7])
        pause(movieLength/length(Traj));
        title(['Particle ID ' num2str(ii)])
    end
    disp(ii);
    %pause;
    close all
end
