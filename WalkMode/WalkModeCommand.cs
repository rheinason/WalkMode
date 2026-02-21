using Rhino;
using Rhino.Commands;

namespace WalkMode
{
    public class WalkModeCommand : Command
    {
        public WalkModeCommand()
        {
            Instance = this;
        }

        public static WalkModeCommand Instance { get; private set; }

        public override string EnglishName => "WalkMode";

        protected override Result RunCommand(RhinoDoc doc, RunMode mode)
        {
            var view = doc.Views.ActiveView;
            if (view == null)
            {
                RhinoApp.WriteLine("No active view.");
                return Result.Failure;
            }

            using (var walker = new WalkGetter(doc, view))
            {
                var walkResult = walker.RunWalk();
                return walkResult == Rhino.Input.GetResult.Cancel
                    ? Result.Cancel
                    : Result.Success;
            }
        }
    }
}
