
%Function applies the kalman filter to sensor data and utilizes the QUEST 
%algorithm to determine the rotation matrix and quaternions representing
%the relative rotations between two reference frames. It also plots the 
%trajectory of the body undergoing the rotation

function [quat, kalmanQuat, refFrameData, kalFlitSensorData] =...
    firstModel(motionData, covMatData, inertialFrameData, version, side, debug)



dataflds = fields(motionData);


for i=1:length(dataflds)
    
    limb = dataflds{i};        
        
        % Filter accelerometer and magnetometer sensors
        [kalFlitSensorData.(limb).g, kalFlitSensorData.(limb).m, kalFlitSensorData.(limb).eps,kalFlitSensorData.(limb).K,kalFlitSensorData.(limb).Ktrace] =...
            kalmanFilter(motionData.(limb).g, motionData.(limb).m, motionData.(limb).w, motionData.(limb).dw, covMatData.(limb),version);
        
        % Rotation estimate with respect to Npose local frame
        
        refFrameData.(limb).g = mean(inertialFrameData.(limb).g')';
        refFrameData.(limb).m = mean(inertialFrameData.(limb).m')';
        
        for j=1:size(motionData.(limb).g,2)
%             j          
            quatIn.(limb).g = motionData.(limb).g(:,j);
            quatIn.(limb).m = motionData.(limb).m(:,j);
            
            kalQuatIn.(limb).g = kalFlitSensorData.(limb).g(:,j);
            kalQuatIn.(limb).m = kalFlitSensorData.(limb).m(:,j);
            
            dataStandDev.v_g =1*100;
            dataStandDev.v_m =1;
            dataStandDev.w_g =1*100;
            dataStandDev.w_m =1;
            
           [kalmanQuat.(fldname).q_opt(:,j),kalmanQuat.(fldname).R{j}, kalmanQuat.(fldname).res_UA{j}] =...
               questalg(kalQuatIn.(fldname),ref.(fldname), dataStandDev);
             

            [quat.(limb).q_opt_q(:,j),quat.(limb).R{j}, quat.(limb).res_UA{j},quat.(limb).q{j}] =...
                questalg(quatIn.(limb),ref.(fldname),dataStandDev);
            
           
                  
        quat.(limb).q_ref = [0 0 0 1]';
        kalmanQuat.(limb).q_ref = [0 0 0 1]';
        
end


end


%plot output
if strcmp(side,'right')
    upfld = 'R_upper_arm';
    lofld = 'R_lower_arm';
else
    upfld = 'L_upper_arm';
    lofld = 'L_lower_arm';
end


if debug
    
    figure
    subplot(3,2,1)
    hold on
    plot(motionData.(lofld).g(1,:),'g')
    plot(kalFlitSensorData.(lofld).g(1,:),'k'),grid

    subplot(3,2,2)
    hold on
    plot(motionData.(lofld).m(1,:),'g')
    plot(kalFlitSensorData.(lofld).m(1,:),'k'),grid

    subplot(3,2,3)
    hold on
    plot(motionData.(lofld).g(2,:),'g')
    plot(kalFlitSensorData.(lofld).g(2,:),'k'),grid

    subplot(3,2,4)
    hold on
    plot(motionData.(lofld).m(2,:),'g')
    plot(kalFlitSensorData.(lofld).m(2,:),'k'),grid

    subplot(3,2,5)
    hold on
    plot(motionData.(lofld).g(3,:),'g')
    plot(kalFlitSensorData.(lofld).g(3,:),'k'),grid

    subplot(3,2,6)
    hold on
    plot(motionData.(lofld).m(3,:),'g')
    plot(kalFlitSensorData.(lofld).m(3,:),'k'),grid



    figure
    for i=1:4
        subplot(4,1,i)
        plot(kalmanQuat.(lofld).q_opt(i,:),'b'),grid
        hold on
        plot(quat.(lofld).q_opt(i,:),'r')
    end



    figure
    subplot(3,2,1)
    hist(kalFlitSensorData.(lofld).eps(1,:))
    title('g')
    subplot(3,2,3)
    hist(kalFlitSensorData.(lofld).eps(2,:))
    subplot(3,2,5)
    hist(kalFlitSensorData.(lofld).eps(3,:))
    subplot(3,2,2)
    title('m')
    hist(kalFlitSensorData.(lofld).eps(4,:))
    subplot(3,2,4)
    hist(kalFlitSensorData.(lofld).eps(5,:))
    subplot(3,2,6)
    hist(kalFlitSensorData.(lofld).eps(6,:))
    
end