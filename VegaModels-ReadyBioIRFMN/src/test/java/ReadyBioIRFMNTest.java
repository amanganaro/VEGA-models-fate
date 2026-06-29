import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.readybio_irfmn.ismReadyBioIRFMN;
import model.ModelExecutionTest;

public class ReadyBioIRFMNTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismReadyBioIRFMN();
    }
}
