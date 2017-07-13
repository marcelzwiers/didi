function [S LWarn] = mywarning(varargin)

% Change the warning state without changing the last warning
%
% SYNTAX
% 
% MYWARNING can be used in the same way as warning, but in addition it can
% capture and restore the warning state together with the last warning. This
% behaviour can be evoked with with the following usage:
%
% [S LWarn] = mywarning(state, 'message_id')
% [..]
% mywarning(S, LWarn)
%
% EXAMPLE
% 
% disp('NORMAL WARNING BEHAVIOUR:')
% warning('This is some previous warning')
% S = warning('Off','test:warningbehaviour');
% warning('test:warningbehaviour','This switched off warning is suppressed but still caught by lastwarn')
% warning(S)
% disp(['Lastwarn: ' lastwarn])
% 
% disp('MYWARNING BEHAVIOUR:')
% warning('This is some previous warning')
% [S LWarn] = mywarning('Off','test:mywarningbehaviour');
% warning('test:mywarningbehaviour','This switched off warning is suppressed and not rethrown by lastwarn')
% mywarning(S, LWarn)
% disp(['Lastwarn: ' lastwarn])
% 
% See also: warning
%
% Marcel, 6-5-2014


if ischar(varargin{1}) && any(strcmpi(varargin{1},{'On','Off','Query'}))
	
	% USAGE: [S LWarn] = mywarning(State, ..)
	% -> Capture the last warning (LWarn) and change the warning state (S)
	[Msg ID] = lastwarn;
	LWarn	 = {Msg ID};
	S		 = warning(varargin{:});
	
elseif isstruct(varargin{1}) && iscell(varargin{2})
	
	% USAGE: mywarning(S, LWarn)
	% -> Restore the warning state (S) and restore the last warning (LWarn)
	warning(varargin{1})					% Restore S
	lastwarn(varargin{2}{:})				% Restore LWarn
	
else
	
	% USAGE: Normal warning usage
	if nargout
		S = warning(varargin{:});
	else
		warning(varargin{:})
	end
	
end
