#ifndef TB_CLIPBOARD_H
#define TB_CLIPBOARD_H 1 

class TbClipboard {
 public:
  TbClipboard ( int nevts ) {m_nevts = nevts};
  ~TbClipboard ();


 private:
  int m_nevts;
