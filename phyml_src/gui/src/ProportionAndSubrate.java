
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;

/**
 * Implements all components for specifying the "Proportion invariable sites"
 * and "One Substitutuon Rate Category" variable.
 * 
 * @author Christoph Knapp
 * 
 */

public class ProportionAndSubrate extends JPanel implements ActionListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private JComboBox sites;
	private CustomTextField parameter;
	private JComboBox subRates;

	/**
	 * Constructor method.
	 */
	public ProportionAndSubrate() {
		String[] choices = new String[] { "fixed", "estimated" };
		sites = new JComboBox(choices);
		sites.addActionListener(this);
		parameter = new CustomTextField("0.00");
		String[] yesNo = new String[] { "yes", "no" };
		subRates = new JComboBox(yesNo);
		subRates.addActionListener(this);
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		JPanel p1 = new JPanel();
		layout.setDimensions(0.45, 0.9);
		add(p1);
		CustomGridLayout lO1 = new CustomGridLayout();
		p1.setLayout(lO1);
		lO1.setDimensions(0.66, 1);
		p1.add(new JLabel("Proportion of Invariable Sites"));
		lO1.setDimensions(0.22, 1);
		p1.add(sites);
		lO1.setDimensions(0.02, 1);
		p1.add(new JPanel());
		lO1.setDimensions(0.1, 1);
		p1.add(parameter);
		JPanel p2 = new JPanel();
		layout.setDimensions(0.04, 0.9);
		add(new JPanel());
		layout.setDimensions(0.49, 0.9);
		add(p2);
		CustomGridLayout lO2 = new CustomGridLayout();
		p2.setLayout(lO2);
		lO2.setDimensions(0.74, 1);
		// p2.add(new JLabel("One Sub. Rate Category"));
		p2.add(new JLabel("Variation of rates across sites"));
		//lO2.setDimensions(0.04, 1);
		//p2.add(new JPanel());
		lO2.setDimensions(0.2, 1);
		p2.add(subRates);
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
	}

	/**
	 * CustomTextField extends JTextField. It is customized that the user is
	 * only able to type in Double values.
	 * 
	 * @author Christoph Knapp
	 * 
	 */
	class CustomTextField extends JTextField {
		/**
		 * default id
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * Constructor Method.
		 * 
		 * @param string
		 *            String : default value displayed in the text field.
		 */
		public CustomTextField(String string) {
			super(string);
			if (!isNumber(string)) {
				this.setText("");
			}
		}

		@Override
		protected Document createDefaultModel() {
			return new doubleOnlyDocument();
		}

		/**
		 * Tests wheter the input String is convertable to a Number.
		 * 
		 * @param str
		 * 
		 * @return
		 */
		private boolean isNumber(String str) {
			try {
				Double.parseDouble(str);
				return true;
			} catch (NumberFormatException e) {
				if (str.equals(".") && !this.getText().contains(".")) {
					return true;
				}
				return false;
			}
		}

		/**
		 * Customizes the default PainDocument to only accept numbers and one
		 * ".". Therefore only Double values can be typed in.
		 * 
		 * @author Christoph Knapp
		 * 
		 */
		class doubleOnlyDocument extends PlainDocument {
			/**
			 * default id
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public void insertString(int offs, String str, AttributeSet a)
					throws BadLocationException {
				if (str == null || !isNumber(str)) {
					return;
				}
				super.insertString(offs, str, a);
			}
		}
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == sites) {
			if (sites.getSelectedItem().toString().equals("fixed")) {
				parameter.setEnabled(true);
			} else {
				parameter.setEnabled(false);
			}
		} else if (e.getSource() == subRates) {
			if (subRates.getSelectedItem().toString().equals("yes")) {
				PhymlPanel.sM.setNumSubCatGammaAverageCompVisible(true);
			} else {
				PhymlPanel.sM.setNumSubCatGammaAverageCompVisible(false);
			}
		}
	}

	/**
	 * Retrieves the value for the command line value for the proportion of
	 * invariable sites.
	 * 
	 * @return String : "e" if estimated or Double value in String format if
	 *         fixed.
	 */
	public String getPropInvarSites() {
		if (sites.getSelectedItem().toString().equals("estimated")) {
			return "e";
		} else {
			return Double.parseDouble(parameter.getText()) + "";
		}
	}
}
