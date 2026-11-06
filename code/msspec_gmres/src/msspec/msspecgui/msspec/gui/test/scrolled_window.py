
import wx

class ScrolledImageComponent(wx.ScrolledWindow):
    
    def __init__(self, parent, id, image_path):
        wx.ScrolledWindow.__init__(self, parent, id)
        
        image = wx.Image(image_path)
        if  not image.IsOk():
            wx.MessageBox("there was an error loading the image")
            return
 
        self.w = image.GetWidth()
        self.h = image.GetHeight()
 
        # init scrolled area size, scrolling speed, etc. */
        self.SetScrollbars(1,1, self.w, self.h, 0, 0)
 
        self.bitmap = wx.Bitmap(image_path)

    def OnDraw(self, dc):
        """
        render the image - in a real app, if your scrolled area
           is somewhat big, you will want to draw only visible parts,
           not everything like below
        """
        dc.DrawBitmap(self.bitmap, 0, 0, False)
 
        #// also check wxScrolledWindow::PrepareDC

# Run the program
if __name__ == "__main__":
    app = wx.App(False)
    wx.InitAllImageHandlers();
    frame = wx.Frame(None)
    sizer = wx.BoxSizer(wx.HORIZONTAL)
    my_image = ScrolledImageComponent(frame, wx.ID_ANY, "/Users/graffy/data/Perso/mms_img1041085029.jpg")
    sizer.Add(my_image, 1, wx.ALL | wx.EXPAND, 120)
    frame.SetSizer(sizer);
    frame.Show()
    app.MainLoop()

"""
class ScrolledImageComponent : public wxScrolledWindow
{
    wxBitmap* bitmap;
    int w,h;
public:
    ScrolledImageComponent(wxWindow* parent, wxWindowID id, wxString image_path) : wxScrolledWindow(parent, id)
    {
        wxImage image(image_path);
        if(!image.IsOk())
        {
            wxMessageBox(wxT("there was an error loading the image"));
            return;
        }
 
        w = image.GetWidth();
        h = image.GetHeight();
 
        /* init scrolled area size, scrolling speed, etc. */
        SetScrollbars(1,1, w, h, 0, 0);
 
        bitmap = new wxBitmap( image );
    }
    ~ScrolledImageComponent()
    {
        delete bitmap;
    }
    void OnDraw(wxDC& dc)
    {
        /* render the image - in a real app, if your scrolled area
           is somewhat big, you will want to draw only visible parts,
           not everything like below */
        dc.DrawBitmap(*bitmap, 0, 0, false);
 
        // also check wxScrolledWindow::PrepareDC
    }
};
 
 
class MyApp: public wxApp
{
    wxFrame *frame;
public:
 
    bool OnInit()
    {
        wxInitAllImageHandlers();
        wxBoxSizer* sizer = new wxBoxSizer(wxHORIZONTAL);
        frame = new wxFrame((wxFrame *)NULL, -1,  wxT("Scrolling an Image"), wxPoint(50,50), wxSize(650,650));
 
        ScrolledImageComponent* my_image = new ScrolledImageComponent(frame, wxID_ANY, wxT("my_image.jpg") );
        sizer->Add(my_image, 1, wxALL | wxEXPAND, 120);
        frame->SetSizer(sizer);
 
        frame->Show();
        return true;
    } 
};
 """
 