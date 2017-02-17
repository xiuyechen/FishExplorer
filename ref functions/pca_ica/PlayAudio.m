classdef PlayAudio < handle
%
% Class for playing and visualizing audio files
%
% Brian Moore
% brimoor@umich.edu
%
% November 10, 2016
%

    %
    % Private constants
    %
    properties (GetAccess = private, Constant = true)
        % Audio
        FPS = 15;                   % Graphics playback, frames/sec
        LATENCY_MS = 500;           % Audioplayer latency, in ms
        WINDOW_MS = 50;             % Spectrogram window width, in ms
        
        % GUI
        GUI_NAME  = 'Audio Player'; % GUI name
        Y_PAD = 0.05;               % Ampitude (y-axis) padding
        LINEWIDTH = 1;              % Status linewidth
        FONT_SIZE = 10;             % Font size
        AXIS_GAP = 40;              % Axis gap, in pixels
        LIST_GAP = 20;              % List gap, in pixels
        TOP_GAP = 35;               % Top gap, in pixels
        LEFT_GAP = 60;              % Left gap, in pixels
        LIST_WIDTH = 150;           % List width, in pixels
        BUTTON_HEIGHT = 20;         % Button height, in pixels
        
        % Colors
        WAVE_COLOR = [93, 147, 191] / 255;  % Waveform color
        STATUS_COLOR = [0, 0, 0] / 255;     % Waveform status color
        PLAY_COLOR = [252, 252, 252] / 255; % Active color
        STOP_COLOR = [204, 51, 51] / 255;   % Stop color
    end
    
    %
    % Public GetAccess properties
    %
    properties (GetAccess = public, SetAccess = private)
        isPlaying = false;          % Playback status
        currIdx = nan;              % Current track index
    end
    
    %
    % Private properties
    %
    properties (Access = public)
        % Data
        audio;                      % Audio data
        
        % Graphics
        fig;                        % Figure handle
        ax1;                        % Waveform axis handle
        ax2;                        % Spectrogram axis handle
        ax3;                        % Track list title axis handle
        listh;                      % Audio list handle
        playh;                      % Play button handle
        ph1;                        % Waveform slider handle
        ph2;                        % Spectrofram slider handle
        xLim;                       % Axes x-limits
        yLim1;                      % Waveform y-limits
        yLim2;                      % Spectrogram y-limits
    end
    
    %
    % Public methods
    %
    methods (Access = public)
        %
        % Constructor
        %
        function this = PlayAudio(audio)
            % Syntax:   player = PlayAudio();
            %           player = PlayAudio(audio);
            
            % Parse inputs
            this.audio = audio;
            
            % Initialize GUI
            this.InitializeGUI();
        end
        
        %
        % Select track
        %
        function Select(this,idx)
            idx = max(1,min(idx,numel(this.audio)));
            this.SetTrack(idx);
        end
        
        %
        % Start playback
        %
        function Start(this)
            if ~this.isPlaying
                resume(this.audio(this.currIdx).player);
            end
            this.UpdateButton();
        end
        
        %
        % Stop playback
        %
        function Stop(this)
            if this.isPlaying
                pause(this.audio(this.currIdx).player);
            end
            this.UpdateButton();
        end
        
        %
        % Toggle playback
        %
        function TogglePlayback(this)
            if this.isPlaying
                this.Stop();
            else
                this.Start();
            end
        end
        
        %
        % Reset playback
        %
        function Reset(this)
            stop(this.audio(this.currIdx).player);
            this.UpdateButton();
            this.UpdateStatusLine();
        end
        
        %
        % Close player
        %
        function Close(this)
            try
                this.Reset();
            catch %#ok
                % Empty
            end
            
            try
                delete(this.fig);
            catch %#ok
                delete(gcf);
            end
        end
    end
    
    %
    % Private methods
    %
    methods (Access = private)
        %
        % Handle list selection
        %
        function HandleListSelection(this)
            value = get(this.listh,'Value');
            this.SetTrack(value);
        end
        
        %
        % Handle keypress
        %
        function HandleKeyPress(this,e)
            % Parse keypress
            keyChar = e.Character;
            if isempty(keyChar)
                return;
            end
            
            % Handle keypress
            switch double(keyChar)
                case 13 % Enter
                    % Toggle playback
                    this.TogglePlayback();
                case {28, 8, 127} % {Left arrow, backspace, delete}
                    % Reset playback
                    this.Reset();
            end
        end
        
        %
        % Set playing status
        %
        function SetPlayStatus(this,isPlaying)
            this.isPlaying = isPlaying;
            this.UpdateButton();
            this.UpdateStatusLine();
        end
        
        %
        % Set track
        %
        function SetTrack(this,idx)
            if idx == this.currIdx
                return;
            elseif ~isnan(this.currIdx)
                this.Reset();
                this.audio(this.currIdx).player = [];
            end
            
            % Update track
            this.currIdx = idx;
            
            % Get data
            y = max(-1,min(this.audio(this.currIdx).y,1));
            Fs = this.audio(this.currIdx).Fs;
            Ny = numel(y);
            t = (1:Ny) / Fs;
            
            % Compute spectrogram
            window = this.WINDOW_MS * (Fs / 1000);
            [S, F, T] = spectrogram(y,window,[],[],Fs,'yaxis'); 
            P = 10 * log10(abs(S));
            
            % Update axis limits
            this.xLim = [t(1), t(Ny)];
            this.yLim1 = [-1, 1];
            this.yLim2 = [min(F), max(F)];
            
            % Update waveform
            axLim1 = [this.xLim, this.yLim1 + this.Y_PAD * [-1, 1]];
            set(this.ph1(1),'XData',t,'YData',y);
            axis(this.ax1,axLim1);
            
            % Update spectrogram
            axLim2 = [this.xLim, this.yLim2];
            set(this.ph1(2),'XData',T, ...
                            'YData',F, ...
                            'ZData',zeros(size(P)), ...
                            'CData',P);
            axis(this.ax2,axLim2);
            
            % Construct audio player
            player = audioplayer(y,Fs);
            player.TimerPeriod = round(1000 / this.FPS) / 1000;
            player.StartFcn = @(s,e)this.SetPlayStatus(true);
            player.StopFcn = @(s,e)this.SetPlayStatus(false);
            player.TimerFcn = @(s,e)this.UpdateStatusLine();
            this.audio(this.currIdx).player = player;
        end
        
        %
        % Update track list
        %
        function UpdateTrackList(this)
            set(this.listh,'String',{this.audio.name});
        end
        
        %
        % Update button
        %
        function UpdateButton(this)
            if this.isPlaying
                set(this.playh,'String','Stop' ,...
                               'BackgroundColor',this.STOP_COLOR, ...
                               'ForegroundColor',[1, 1, 1]);
            else
                set(this.playh,'String','Play', ...
                               'BackgroundColor',this.PLAY_COLOR, ...
                               'ForegroundColor',[0, 0, 0]);
            end
        end
        
        %
        % Update status line
        %
        function UpdateStatusLine(this)
            % Get current time
            if this.isPlaying
                n = this.audio(this.currIdx).player.CurrentSample;
                Fs = this.audio(this.currIdx).Fs;
            else
                n = -1;
                Fs = 1;
            end
            t = (n / Fs) - (this.LATENCY_MS / 1000);
            
            % Update waveform status line
            x1 = [t, t];
            y1 = this.yLim1;
            set(this.ph2(1),'XData',x1,'YData',y1);
            
            % Update spectrogram status line
            x2 = [t, t];
            y2 = this.yLim2;
            set(this.ph2(2),'XData',x2,'YData',y2);
        end
        
        %
        % Initialize GUI
        %
        function InitializeGUI(this)
            % Setup figure
            this.fig = figure( ...
                'Name',this.GUI_NAME, ...
                'MenuBar','none', ...
                'NumberTitle','off', ...
                'DockControl','off', ...
                'KeyPressFcn',@(s,e)this.HandleKeyPress(e), ...
                'ResizeFcn',@(s,e)this.ResizeGUI(), ...
                'CloseRequestFcn',@(s,e)this.Close(), ...
                'Visible','off');
            
            % File menu
            filem = uimenu(this.fig,'Label','File');
            uimenu(filem,'Label','Close', ...
                         'Callback',@(s,e)this.Close(), ...
                         'Accelerator','W');
            
            % Play menu
            playm = uimenu(this.fig,'Label','Play');
            uimenu(playm,'Label','Start', ...
                         'Callback',@(s,e)this.Start(), ...
                         'Accelerator','1');
            uimenu(playm,'Label','Stop', ...
                         'Callback',@(s,e)this.Stop(), ...
                         'Accelerator','2');
            uimenu(playm,'Label','Reset', ...
                         'Callback',@(s,e)this.Reset(), ...
                         'Accelerator','3', ...
                         'Separator','on');
            
            % Track list
            this.ax3 = axes('Units','pixels');
            title(this.ax3,'Tracks');
            axis(this.ax3,'off');
            this.listh = uicontrol( ...
                'Parent',this.fig, ...
                'Style','List', ...
                'Value',1, ...
                'Max',1, ...
                'Enable','on', ...
                'Units','pixels', ...
                'HorizontalAlignment','left', ...
                'FontSize',this.FONT_SIZE, ...
                'Callback',@(s,e)this.HandleListSelection(), ...
                'KeyPressFcn',@(s,e)this.HandleKeyPress(e));
            
            % Play button
            this.playh = uicontrol( ...
                      'Parent',this.fig, ...
                      'Style','pushbutton', ...
                      'String','', ...
                      'Units','pixels', ...
                      'FontSize',this.FONT_SIZE, ...
                      'HorizontalAlignment','center', ...
                      'Callback',@(s,e)this.TogglePlayback());
            
            % Waveform
            this.ax1 = axes('Units','pixels');
            this.ph1(1) = plot( ...
                this.ax1,nan,nan,'-','Color',this.WAVE_COLOR);
            ylabel(this.ax1,'Amplitude','FontSize',this.FONT_SIZE);
            title(this.ax1,'Waveform','FontSize',this.FONT_SIZE);
            set(this.ax1,'FontSize',this.FONT_SIZE);
            axis(this.ax1,'manual');
            hold(this.ax1,'on');
            
            % Spectrogram
            this.ax2 = axes('Units','pixels');
            this.ph1(2) = pcolor(this.ax2,nan(2));
            xlabel(this.ax2,'Time (s)','FontSize',this.FONT_SIZE);
            ylabel(this.ax2,'Frequency (Hz)','FontSize',this.FONT_SIZE);
            title(this.ax2,'Spectrogram','FontSize',this.FONT_SIZE);
            set(this.ax2,'FontSize',this.FONT_SIZE);
            shading(this.ax2,'flat');
            hold(this.ax2,'on');
            
            % Waveform status line
            this.ph2(1) = plot( ...
                nan,nan,'-', ...
                'Color',this.STATUS_COLOR, ...
                'Linewidth',this.LINEWIDTH, ...
                'Parent',this.ax1);
            
            % Spectrogram waveform status line
            this.ph2(2) = plot( ...
                nan,nan,'-', ...
                'Color',this.STATUS_COLOR, ...
                'Linewidth',this.LINEWIDTH, ...
                'Parent',this.ax2);
            
            % Initialize GUI
            this.UpdateButton();
            this.UpdateTrackList();
            this.SetTrack(1);
            this.ResizeGUI();
            set(this.fig,'Visible','on');
        end
        
        %
        % Resize GUI
        %
        function ResizeGUI(this)
            % Get constants
            gapa = this.AXIS_GAP;
            gapb = this.LIST_GAP;
            gapt = this.TOP_GAP;
            gapl = this.LEFT_GAP;
            w2 = this.LIST_WIDTH;
            bh = this.BUTTON_HEIGHT;
            
            % Get figure dimensions
            pos = get(this.fig,'Position');
            xt = pos(3);
            yt = pos(4);
            
            % Compute GUI element positions
            w1 = xt - w2 - gapl - 2 * gapb;
            h1 = 0.5 * (yt - 2 * gapa - gapt);
            h2 = yt - gapa - gapb - gapt - bh;
            ax1Pos = [gapl, (2 * gapa + h1), w1, h1];
            ax2Pos = [gapl, gapa, w1, h1];
            listPos = [(gapl + gapb + w1), (gapa + gapb + bh), w2, h2];
            playPos = [(gapl + gapb + w1), gapa, w2, bh];
            
            % Update positions
            set(this.ax1,'Position',ax1Pos);
            set(this.ax2,'Position',ax2Pos);
            set(this.ax3,'Position',listPos);
            set(this.listh,'Position',listPos);
            set(this.playh,'Position',playPos);
        end
    end
end
