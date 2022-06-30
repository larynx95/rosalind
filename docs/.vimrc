call plug#begin()

Plug 'tpope/vim-surround'
Plug 'vim-airline/vim-airline'
Plug 'kien/ctrlp.vim'
Plug 'itchyny/landscape.vim'
Plug 'scrooloose/nerdtree'
Plug 'scrooloose/syntastic'

call plug#end()

set encoding=UTF-8
set nu
set ruler
set vb
set nobackup
set nowritebackup
set noswapfile
set ts=4
set shiftwidth=4
set linespace=3
colorscheme landscape
highlight Normal guibg=black guifg=white
set background=dark

"Mode Settings
let &t_SI.="\e[5 q" "SI = INSERT mode
let &t_SR.="\e[4 q" "SR = REPLACE mode
let &t_EI.="\e[1 q" "EI = NORMAL mode (ELSE)
"Cursor settings:
"  1 -> blinking block
"  2 -> solid block
"  3 -> blinking underscore
"  4 -> solid underscore
"  5 -> blinking vertical bar
"  6 -> solid vertical bar

" In insert or command mode, move normally by using Ctrl
inoremap <C-h> <Left>
inoremap <C-j> <Down>
inoremap <C-k> <Up>
inoremap <C-l> <Right>
cnoremap <C-h> <Left>
cnoremap <C-j> <Down>
cnoremap <C-k> <Up>
cnoremap <C-l> <Right>
nnoremap <C-n> :NERDTreeFocus<CR>

" syntastic
set statusline+=%#warningmsg#
set statusline+=%{SyntasticStatuslineFlag()}
set statusline+=%*
let g:syntastic_always_populate_loc_list = 1
let g:syntastic_auto_loc_list = 1
let g:syntastic_check_on_open = 1
let g:syntastic_check_on_wq = 0

" autocmd
autocmd FileType python map <buffer> <F5> :w<CR> :cd %:p:h<CR> :exec '!python3' shellescape(@%, 1)<CR>
autocmd FileType go map <buffer> <F5> :w<CR> :cd %:p:h<CR> :exec '!go run' shellescape(@%, 1)<CR>
autocmd FileType javascript map <buffer> <F5> :w<CR> :cd %:p:h<CR> :exec '!node' shellescape(@%, 1)<CR>
autocmd filetype cpp nnoremap <f5> :w<CR> :cd %:p:h<CR> :exec '!g++ -std=c++2a -Wall % -o %:r && ./%:r' <CR> :exec '!rm ./%:r'
autocmd filetype rust nnoremap <f5> :w<CR> :cd %:p:h<CR> :exec '!cargo run' shellescape(@%, 1)<CR> :exec '!cargo clean'
