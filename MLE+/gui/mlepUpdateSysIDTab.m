function [handles] = mlepUpdateSysIDTab(handles)

%% EXISTENCE
% Variable Input
if isfield(handles.DATA, 'variableInput') && ~isempty(handles.DATA.variableInput) && ~isfield(handles.DATA, 'SystemID_InputListbox')
    handles.DATA.SystemID_InputListbox = handles.DATA.variableInput(:,5);
    handles.DATA.SystemID_InputCommentEdit = handles.DATA.variableInput(:,4);
end

% Variable Output
if isfield(handles.DATA, 'variableOutput') && ~isempty(handles.DATA.variableOutput) && ~isfield(handles.DATA, 'SystemID_OutputListbox')
    handles.DATA.SystemID_OutputListbox = handles.DATA.variableOutput(:,5);
    handles.DATA.SystemID_OutputCommentEdit = handles.DATA.variableOutput(:,4);
end

% Input Listbox
if ~isfield(handles.DATA, 'SystemID_InputListbox')
    handles.DATA.SystemID_InputListbox = [];
    handles.DATA.SystemID_InputCommentEdit = [];
    set(handles.SystemID_InputListbox, 'String', handles.DATA.SystemID_InputListbox);
    set(handles.SystemID_InputCommentEdit, 'String', handles.DATA.SystemID_InputCommentEdit);
else
    set(handles.SystemID_InputListbox, 'String', handles.DATA.SystemID_InputListbox);
    set(handles.SystemID_InputCommentEdit, 'String', handles.DATA.SystemID_InputCommentEdit);
end

% Output Listbox
if ~isfield(handles.DATA, 'SystemID_OutputListbox')
    handles.DATA.SystemID_OutputListbox = [];
    handles.DATA.SystemID_OutputCommentEdit = [];
    set(handles.SystemID_OutputListbox, 'String', handles.DATA.SystemID_OutputListbox);
    set(handles.SystemID_OutputCommentEdit, 'String', handles.DATA.SystemID_OutputCommentEdit);
else
    set(handles.SystemID_OutputListbox, 'String', handles.DATA.SystemID_OutputListbox);
    set(handles.SystemID_OutputCommentEdit, 'String', handles.DATA.SystemID_OutputCommentEdit);
end

% Control File
if ~isfield(handles.DATA,'ControlFileName')
    handles.DATA.ControlFileName = [];
    set(handles.SystemID_LoadControlFileEdit, 'String', 'Control File');
    set(handles.SystemID_CreateControlFileEdit, 'String', 'ControlFile.m');
    set(handles.SystemID_LoadControlFileEdit, 'Background', 'white');
    set(handles.SystemID_CreateControlFileEdit, 'Background', 'white');
else
    if isfield(handles.DATA,'ControlFileCreated')
        if handles.DATA.ControlFileCreated == 1
            set(handles.SystemID_CreateControlFileEdit, 'String', handles.DATA.ControlFileName);
            set(handles.SystemID_CreateControlFileEdit, 'Background', 'c');
            set(handles.SystemID_LoadControlFileEdit, 'Background', 'white');
        else
            set(handles.SystemID_LoadControlFileEdit, 'String', handles.DATA.ControlFileName);
            set(handles.SystemID_LoadControlFileEdit, 'Background', 'c');
            set(handles.SystemID_CreateControlFileEdit, 'Background', 'white');
        end
    end
end


% %% EXISTENCE
% % Input Listbox
% if ~isfield(handles.DATA, 'SystemID_InputListbox')
%     handles.DATA.SystemID_InputListbox = [];
%     handles.DATA.SystemID_InputCommentEdit = [];
%     set(handles.SystemID_InputListbox, 'String', handles.DATA.SystemID_InputListbox);
%     set(handles.SystemID_InputCommentEdit, 'String', handles.DATA.SystemID_InputCommentEdit);
% end
% 
% % Output Listbox
% if ~isfield(handles.DATA, 'SystemID_OutputListbox')
%     handles.DATA.SystemID_OutputListbox = [];
%     handles.DATA.SystemID_OutputCommentEdit = [];
%     set(handles.SystemID_OutputListbox, 'String', handles.DATA.SystemID_OutputListbox);
%     set(handles.SystemID_OutputCommentEdit, 'String', handles.DATA.SystemID_OutputCommentEdit);
% end
% 
% % Control File
% if ~isfield(handles.DATA,'SystemID_ControlFileName')
% 	handles.DATA.SystemID_ControlFileName = [];
%     set(handles.SystemID_LoadControlFileEdit, 'String', 'Control File');
%     set(handles.SystemID_CreateControlFileEdit, 'String', 'ControlFile.m');
%     set(handles.SystemID_LoadControlFileEdit, 'Background', 'white');
%     set(handles.SystemID_CreateControlFileEdit, 'Background', 'white');
% end
% 
% % Workspace Vars
% if ~isfield(handles.DATA,'workVars')
% 	handles.DATA.workVars = [];
% end


%% CHECK SIZE - Inputs/Input ID Listbox - Comments
% if size(handles.DATA.listInputExt2,1)
%     % INput Variables
%     set(handles.sysIDInComment, 'String',handles.DATA.listInputExt2{handles.DATA.sysIDListInputIndex,6}, 'Enable', 'on');
%     set(handles.sysIDInListbox, 'Value', handles.DATA.sysIDListInputIndex);
%     set(handles.sysIDInListbox, 'String', handles.DATA.listInputExt2(:,7), 'Enable', 'on');
%     set(handles.sysIDInListboxID, 'Value', handles.DATA.sysIDListInputIndex);
%     set(handles.sysIDInListboxID, 'String', handles.DATA.listInputExt2(:,7), 'Enable', 'on');
% else
%     % Inactive
%     set(handles.sysIDInComment, 'String', [], 'Enable', 'off');
%     set(handles.sysIDInListbox, 'String', [], 'Enable', 'off');
%     set(handles.sysIDInListboxID, 'String', [], 'Enable', 'off');
% end

%% CHECK SIZE - Outputs Listbox - Comments
% if size(handles.DATA.listOutputExt2,1)
%     % OUTput Variables
%     set(handles.sysIDOutComment, 'String',handles.DATA.listOutputExt2{handles.DATA.sysIDListOutputIndex,5}, 'Enable', 'on');
%     set(handles.sysIDOutListbox, 'Value', handles.DATA.sysIDListOutputIndex);
%     set(handles.sysIDOutListbox, 'String', handles.DATA.listOutputExt2(:,6), 'Enable', 'on');
% else
%     % Inactive
%     set(handles.sysIDOutComment, 'String', [], 'Enable', 'off');
%     set(handles.sysIDOutListbox, 'String', [], 'Enable', 'off');
% end

%% Set Variables Listbox
% if ~isfield(handles.DATA,'sysIDvarsIndex')
%     handles.DATA.sysIDvarsIndex = 1;
% end
% set(handles.sysIDvarsListbox, 'Value', handles.DATA.sysIDvarsIndex);
% 
% if isfield(handles.DATA,'sysIDvarsText')
%     set(handles.sysIDvarsListbox, 'String', handles.DATA.sysIDvarsText, 'Enable', 'on');
% else
%     set(handles.sysIDvarsListbox, 'Enable', 'off');
% end

%% Set Input IDDATA
% Set index of input listbox 
% if ~isfield(handles.DATA,'sysIDInputTextSelectedIndex')
%     handles.DATA.sysIDInputTextSelectedIndex = 1;
% end
% set(handles.sysIDInputListboxSelected, 'Value', handles.DATA.sysIDInputTextSelectedIndex);
% 
% % Check for Input selected
% if isfield(handles.DATA,'sysIDInputTextSelected')
%     set(handles.sysIDInputListboxSelected, 'String', handles.DATA.sysIDInputTextSelected, 'Enable', 'on');
% else
%     set(handles.sysIDInputListboxSelected, 'Enable', 'off');
% end

%% Set Output IDDATA
% Set index of output listbox 
% if ~isfield(handles.DATA,'sysIDOutputTextSelectedIndex')
%     handles.DATA.sysIDOutputTextSelectedIndex = 1;
% end
% set(handles.sysIDOutputListboxSelected, 'Value', handles.DATA.sysIDOutputTextSelectedIndex);
% 
% % Check for Output selected
% if isfield(handles.DATA,'sysIDOutputTextSelected')
%     set(handles.sysIDOutputListboxSelected, 'String', handles.DATA.sysIDOutputTextSelected, 'Enable', 'on');
% else
%     set(handles.sysIDOutputListboxSelected, 'Enable', 'off');
% end

%% SET SYS PARAM
[handles] = mlepUpdateSysIDParam(handles);
end