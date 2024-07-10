function [Img, Tracks, Tracks_Info] = read_video(A,rootfldr ,Img_name, No_frms, No_tracks)

    Imfrmt = '.jpg';
    if No_frms < 100
        str1 = '_01'; 
        str2 = '%02d'; 
%         str1 = '_324'; 
%         str2 = '%03d'; 
    else
        str1 = '_001'; 
%         str1 = '454'; 
        str2 = '%03d';
    end;

    Im_sz = size(imread([rootfldr Img_name '/' Img_name str1 Imfrmt])); 
    Img = zeros(Im_sz(1), Im_sz(2), 3*No_frms); 
    for i=1:No_frms
        Img(:, :, (i-1)*3+1:i*3) = double(imread([rootfldr Img_name '/' Img_name '_' num2str(i, str2) Imfrmt])); 
%         Img(:, :, (i-1)*3+1:i*3) = double(imread([rootfldr Img_name '/' Img_name '_' num2str(i+324-1, str2) Imfrmt])); 
%         Img(:, :, (i-1)*3+1:i*3) = double(imread([rootfldr Img_name '/' Img_name num2str(i+454-1, str2) Imfrmt])); 
    end
    
    count = 3;
    Tracks = zeros(No_tracks, 2, No_frms); 
    Tracks_Info = zeros(6, No_tracks); 

    for i=1:No_tracks

        Label = A(count); 
        count = count+1; 
        Length = A(count); 

        no_sz = 3*Length; 
        Points = reshape(A(count+1:count+no_sz), [3, Length]); 

        Index = Points(3, :); 
        Tracks(i, :, Index(1)+1:Index(end)+1) = Points(1:2, :); 

        Tracks_Info(1:2, i) =    mean(Points(1:2, :), 2);  %mean of each point in frames
        Tracks_Info(3, i) = mean(Index);        
        Tracks_Info(4, i) = Index(1)+1;         %start of Track
        Tracks_Info(5, i) = Index(end)+1;       %end of track
        Tracks_Info(6, i) = Label+1;            %Ground truth label

        count = count + no_sz + 1; 

    end

end