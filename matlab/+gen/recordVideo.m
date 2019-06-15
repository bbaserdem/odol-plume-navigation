function [SPH, ERR] = recordVideo(SPH, NAME, TIME, FRAME)
%RECORDVIDEO Create a video of the 
%   Create a video file detailing the fluid simulation

if ~exist('SPH','var')
    SPH = init.lidDrivenCavity(50,10000);
end
if ~exist('NAME','var')
    NAME = ['data/video_', datestr(now,'dd-mm-yyyy_HH:MM:SS')];
end
if ~exist('TIME','var')
    TIME = 50000;
end
if ~exist('FRAME','var')
    FRAME = 600;
end
PERIOD = max(1, floor(TIME/(FRAME-1)));
VIDEO = VideoWriter(NAME);

% Open video object
open(VIDEO);

% Draw the setup, and get figure handle
SPH.showParticles;
IMG = SPH.figPar;
figure(IMG);

% Write this frame to the video
VIDEO.writeVideo(getframe(IMG));

% Make sure to finish on errors
try
    for t = 1:TIME
        % Do one step of the simulation
        fprintf('Simulation step, %d\n', t);
        SPH.doTimeStep;
        % Will write this frame if requested
        if mod(t, PERIOD) == 0
            fprintf('Writing frame...\n', t);
            % Draw and update figure, and grab handle
            SPH.showParticles;
            drawnow;
            IMG = SPH.figPar;
            % Write this handle to the video
            VIDEO.writeVideo(getframe(IMG));
        end
    end
    % Close the file
    close(VIDEO);
catch ERR
    % Close the file
    close(VIDEO);
    % Rethrow error
    warning(ERR.message);
end

% END OF FUNCTION

end

