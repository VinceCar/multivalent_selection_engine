% prep_data_for_interpolation_in_matlab
particle_status(:,1)=1; % this is k On
if Valency==3
    load(Source_kD_file)
    a0 = 0.4; % in nm
    Xq=zeros(Number_of_particles,5);
    Xq(:,1:3)=(particle(:,:,3)/a0).^(1/0.6)+1; % conver in peg monomers
% Xq(:,1:3)=(particle(:,:,3)/a0); % conver in peg monomers
    Xq(:,4:5)=particle(:,2:3,2)./particle(:,1,2);
%     Xq(:,4:5)=[sqrt(sum(Xq(:,4:5).^2,2)),Xq(:,5)./Xq(:,4)];
%     Xq(:,4)=Xq(:,4).*cos(Xq(:,5));
    Xq(:,5)=log(atan(Xq(:,5) ./ Xq(:,4)) / pi);
    Xq(:,4)=log(Xq(:,4));

    particle_status(:,1)=1;
    y_1=interpn(X_g_1,X_g_2,X_g_3,X_g_4,X_g_5,y_1_g,Xq(:,1),Xq(:,2),Xq(:,3),Xq(:,4),Xq(:,5),'linear');
    y_2=interpn(X_g_1,X_g_2,X_g_3,X_g_4,X_g_5,y_2_g,Xq(:,1),Xq(:,2),Xq(:,3),Xq(:,4),Xq(:,5),'linear');
    y_3=interpn(X_g_1,X_g_2,X_g_3,X_g_4,X_g_5,y_3_g,Xq(:,1),Xq(:,2),Xq(:,3),Xq(:,4),Xq(:,5),'linear');
    b=0.68;
    omegaLU=0.03;
    particle_status(:,2)=[y_1./(omegaLU.*b./a0.*particle(:,1,2))+y_2./((omegaLU.*b./a0.*particle(:,1,2)).^2)+y_3./((omegaLU.*b./a0.*particle(:,1,2)).^3)].^(-1);
else
    particle_status(:,2)=particle(:,1,2); % this is kD
end
particle_status(:,3)=0;
