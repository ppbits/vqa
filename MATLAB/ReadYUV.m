% (c) Mani Malek Esmaeili
% This function reads a YUV file assuming that the file into RGB matrices
% of the same size (as frame_size)

%function [Y U V]=Read_YUV(FilePath,frame_size)
function video = ReadYUV(FilePath, frame_size, color2gray)
if nargin < 3
    color2gray = true;
end

fid=fopen(FilePath,'r');
if fid<0
%    Y=[];
%    U=[];
%    V=[];
	video = [];
else
    fseek(fid,0,'eof');
    len=ftell(fid)/(frame_size(1)*frame_size(2)*1.5);
    fseek(fid,0,'bof');
%     Y=uint8(zeros(frame_size(1),frame_size(2),len));
%     U=Y;
%     V=Y;
    frame=uint8(zeros(frame_size(1),frame_size(2),3));
    
    if color2gray
        video = uint8(zeros(frame_size(1), frame_size(2), len));
    else
        video = uint8(zeros(frame_size(1), frame_size(2), 3, len));
    end
     
	for k=1:len
        temp=fread(fid,[frame_size(2) frame_size(1)],'uint8');
        %         Y(:,:,k)=temp';
        frame(:,:,1)=temp';
        temp=fread(fid,[frame_size(2)/2 frame_size(1)/2],'uint8');
        %         U(1:2:end,1:2:end,k)=temp';
        frame(1:2:end,1:2:end,2)=temp';
        frame(2:2:end,1:2:end,2)=temp';
        frame(1:2:end,2:2:end,2)=temp';
        frame(2:2:end,2:2:end,2)=temp';
        
        temp=fread(fid,[frame_size(2)/2 frame_size(1)/2],'uint8');
        %         V(1:2:end,1:2:end,k)=temp';
        frame(1:2:end,1:2:end,3)=temp';
        frame(2:2:end,1:2:end,3)=temp';
        frame(1:2:end,2:2:end,3)=temp';
        frame(2:2:end,2:2:end,3)=temp';
        rgb=ycbcr2rgb(frame);% convert from RGB to gray scale
        if color2gray
            video(:,:,k) = rgb2gray(rgb);
        else
            video(:,:,:,k) = rgb;
        end
    end
end
