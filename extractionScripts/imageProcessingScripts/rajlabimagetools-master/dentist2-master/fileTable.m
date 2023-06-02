classdef fileTable < handle
    
    properties (Access = public)
        
        files
        channels % should be protected; use get.channels
        
    end
    
    methods
        
        function p = fileTable(varargin)
            if nargin == 0
                fprintf('New Table\n');
                p.files = cell2table(cell(0,7), 'VariableNames', {'fileID', 'fileName', 'channel', 'top', 'left', 'height', 'width'});
            else
                fprintf('Loading Table\n');
                p.files = readtable(varargin{1},'TextType','string');
            end
        end
        
        function p = loadFilesLegacy(p)
            % implement row snake
            % assume overlap of 100px
            nRows = 5;
            nCols = 5;
            overlap = 100; % pixel overlap
            
            fileName = [];
            channel = [];
            top = [];
            left = [];
            height = [];
            width = [];
            
            function addToTable(filePattern)
                dd = dir(filePattern);
                for i = 1:numel(dd)
                    im = imread(dd(i).name);
                    fileName = [fileName string(dd(i).name)];
                    height = [height size(im,1)];
                    width = [width size(im,2)];
                    channel = [channel string(dd(1).name(1:end-7))];
                end
                
                % insert snaking code here
                ht = size(im,1); % height
                wd = size(im,2); % width
                for i = 0:nRows-1
                    top = [top i*(ht-overlap)*ones(1,nCols)+1];
                    if (mod(i+1,2)==1) % odd
                        left = [left (0:nCols-1)*(wd-overlap)+1];
                    else
                        left = [left (nCols-1:-1:0)*(wd-overlap)+1];
                    end
                end
            end
            
            addToTable('dapi*');
            addToTable('tmr*');
            addToTable('alexa*');
            addToTable('trans*');
            
            p.files = table((1:length(fileName))',fileName',channel', top', left', height', width','VariableNames', {'fileID', 'fileName', 'channel', 'top', 'left', 'height', 'width'});
                        
        end
        
        function channels = get.channels(p)
            channels = unique(p.files.channel);
        end
        
        
        function [paddedIm, outRect] = getPaddedImage(p,fileID,padding)
            channel = p.files.channel(fileID);
            file = p.files(p.files.fileID==fileID,:);
            %rect = file.
            rect = [file.top-padding, file.left-padding, file.height+2*padding, file.width+2*padding];
            [paddedIm, outRect] = p.getImageFromRect(rect,channel,"largest");
            
            %[cropIm, outRect] = p.getImageFromRect(
%             files = p.files(p.files.channel == channel); % keep just this channel
            
        end
        
        function [cropIm, outRect] = getImageFromRect(p,rect,channel,option)
            % rect is [top, left, height, width]
            % top and left will be included, and final image is of size
            % height, width
            % option can be 'largest' to keep largest tile or
            % 'acquisition' to favor the first acquired.
            % eg: cropIm = ft.getImageFromRect([1 1 1022 1024],'dapi','largest');
%             i = 1;
            files = p.files(p.files.channel == string(channel),:);
            B = [files.top files.left files.width files.height];
            
            ri = rectint(rect,B);
            idx = find(ri);
            if string(option) == "largest" % default is "acquisition"
                [~,sortidx] = sort(ri(idx),'descend');
                idx = idx(sortidx);
            end
            %idx = find(rectint(rect,B)); % These are the rectangles that overlap the requested rectangle
            
            
            top = min(files.top(idx));
            left = min(files.left(idx));
            bottom = max(files.top(idx) + files.height(idx));
            right = max(files.left(idx) + files.width(idx));
            
            fullIm = zeros(bottom-top+1,right-left+1,'uint16');
            
            for i = flip(idx) % last to first or smallest to biggest
                im = imread(files.fileName(i));
                topIdx = files.top(i)-top+1;
                leftIdx = files.left(i)-left+1;
                fullIm(topIdx:topIdx+files.height(i)-1,leftIdx:leftIdx+files.width(i)-1) = im;
            end
            
            if rect(1) < top
                rect(3) = rect(3)-(top-rect(1));
                rect(1) = top;
            end
            if rect(2) < left
                rect(4) = rect(4)-(left-rect(2));
                rect(2) = left;
            end
            if rect(3) + rect(1) - top > bottom-top+1
                rect(3) = bottom-top+1-rect(1)+top;
            end
            if rect(4) + rect(2) - left > right-left+1
                rect(4) = right-left+1-rect(2)+left;
            end
            
            cropIm = fullIm(rect(1)-top+1:rect(1)-top+rect(3),rect(2)-left+1:rect(2)-left+rect(4));
            
            outRect = rect;
        end
        
        
    end
    
end 